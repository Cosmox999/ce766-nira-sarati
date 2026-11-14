"""Finish the Sarati drought and flow-intermittency analysis."""
# Author: Ganesh Sonawane
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/sarati_source"
OUT = ROOT / "analysis/final_results"
START, END = "1981-06-01", "2011-05-31"
MONTHS = [6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5]
LABELS = ["Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May"]
YEARS = np.arange(1981, 2011)
SEASONS = {"Jun-Sep": [6, 7, 8, 9], "Oct-Nov": [10, 11],
           "Dec-Feb": [12, 1, 2], "Mar-May": [3, 4, 5]}


def read_daily(raw=RAW):
    columns = ["year", "month", "day"]
    p = pd.read_csv(raw / "catchment_mean_forcings/04049.csv")
    q = pd.read_csv(raw / "streamflow_timeseries/streamflow_observed.csv")
    expected = pd.date_range("1980-01-01", "2020-12-31")
    for frame in (p, q):
        frame.index = pd.to_datetime(frame[columns])
        if not frame.index.equals(expected):
            raise ValueError("The input dates are duplicated, missing or out of order.")
    daily = pd.DataFrame({"P_mm": p["prcp(mm/day)"], "Q_m3s": q["4049"]})
    if (daily < 0).any().any() or np.isinf(daily.to_numpy()).any():
        raise ValueError("Check negative or infinite input values.")
    study = daily.loc[START:END]
    if len(study) != 10957 or study.isna().any().any():
        raise ValueError("The stated study window must have 10,957 complete daily pairs.")
    if daily.P_mm.isna().any():
        raise ValueError("Rainfall is incomplete in the SPI reference record.")
    return daily


def monthly_data(daily):
    if not np.allclose(daily.P_mm * 100, np.round(daily.P_mm * 100), rtol=0, atol=1e-9):
        raise ValueError("Update the aggregation precision for rainfall finer than 0.01 mm.")
    monthly = daily.resample("MS").agg(
        P_mm=("P_mm", "sum"), P_n=("P_mm", "count"),
        Q_n=("Q_m3s", "count"), zero_n=("Q_m3s", lambda x: x.eq(0).sum()),
        nearzero_n=("Q_m3s", lambda x: x.le(0.1).sum()))
    monthly.loc[monthly.P_n.ne(monthly.index.days_in_month), "P_mm"] = np.nan
    # Preserve ties at the 0.01 mm precision stored in the source CSV.
    monthly["P_mm"] = monthly.P_mm.round(2)
    monthly["zero_fraction"] = monthly.zero_n / monthly.Q_n.replace(0, np.nan)
    monthly["nearzero_fraction"] = monthly.nearzero_n / monthly.Q_n.replace(0, np.nan)
    return monthly


def gamma_parameters(values, method="mle"):
    positive = np.asarray(values)[np.asarray(values) > 0]
    if len(positive) < 3 or np.ptp(positive) == 0:
        raise ValueError("Too few distinct positive totals for a gamma fit.")
    if method == "mle":
        shape, _, scale = stats.gamma.fit(positive, floc=0)
    elif method == "moments":
        mean, variance = positive.mean(), positive.var(ddof=0)
        shape, scale = mean * mean / variance, variance / mean
    else:
        raise ValueError(method)
    return float(shape), float(scale)


def gamma_spi(values, shape, scale, p_zero):
    probability = p_zero + (1 - p_zero) * stats.gamma.cdf(values, shape, scale=scale)
    return stats.norm.ppf(np.clip(probability, 1e-10, 1 - 1e-10))


def calculate_spi(monthly, reference=(1981, 2020), method="mle"):
    output = monthly.copy()
    fits = []
    for k in (3, 6):
        total = monthly.P_mm.rolling(k, min_periods=k).sum().round(2)
        output[f"R{k}_mm"] = total
        output[f"SPI{k}"] = np.nan
        for month in range(1, 13):
            ref = total[(total.index.year >= reference[0]) &
                        (total.index.year <= reference[1]) & (total.index.month == month)].dropna()
            if len(ref) != reference[1] - reference[0] + 1:
                raise ValueError("The SPI reference sample is incomplete.")
            shape, scale = gamma_parameters(ref.to_numpy(), method)
            p_zero = float(ref.eq(0).mean())
            use = total[(total.index.month == month) & total.notna()]
            output.loc[use.index, f"SPI{k}"] = gamma_spi(use, shape, scale, p_zero)
            fits.append(dict(spi_scale=k, ending_month=month, n_reference=len(ref),
                             positive_n=int(ref.gt(0).sum()), zero_probability=p_zero,
                             shape=shape, scale_mm=scale, reference_start=reference[0],
                             reference_end=reference[1], fitting_method=method))
    return output, pd.DataFrame(fits)


def drought_events(spi, study_start="1981-06-01", study_end="2011-05-01"):
    """Negative sequences reaching -1, with boundary truncation recorded."""
    study = spi.loc[study_start:study_end]
    rows, run = [], []

    def finish():
        if not run:
            return
        values = study.loc[run]
        if values.min() > -1:
            return
        first, last = run[0], run[-1]
        previous, following = first - pd.offsets.MonthBegin(1), last + pd.offsets.MonthBegin(1)
        rows.append(dict(start=first.strftime("%Y-%m"), end=last.strftime("%Y-%m"),
                         duration_months=len(run), minimum_spi=float(values.min()),
                         peak_month=values.idxmin().strftime("%Y-%m"),
                         magnitude=float(-values.sum()),
                         left_censored=bool(first == study.index[0] and
                                            (previous not in spi.index or spi.loc[previous] < 0)),
                         right_censored=bool(last == study.index[-1] and
                                             (following not in spi.index or spi.loc[following] < 0))))

    for date, value in study.items():
        if pd.notna(value) and value < 0:
            run.append(date)
        else:
            finish()
            run = []
    finish()
    return pd.DataFrame(rows, columns=["start", "end", "duration_months", "minimum_spi",
                                      "peak_month", "magnitude", "left_censored", "right_censored"])


def block_indices(rng, n, repetitions, length=1):
    """Draw whole years, or circular runs of two or three adjacent years."""
    starts = rng.integers(0, n, size=(repetitions, int(np.ceil(n / length))))
    return ((starts[:, :, None] + np.arange(length)) % n).reshape(repetitions, -1)[:, :n]


def row_spearman(x, y):
    """Re-rank each resample, including its repeated years and tied values."""
    xr = stats.rankdata(x, axis=1, method="average")
    yr = stats.rankdata(y, axis=1, method="average")
    xr -= xr.mean(axis=1, keepdims=True)
    yr -= yr.mean(axis=1, keepdims=True)
    denominator = np.sqrt(np.sum(xr * xr, axis=1) * np.sum(yr * yr, axis=1))
    return np.divide(np.sum(xr * yr, axis=1), denominator,
                     out=np.full(len(xr), np.nan), where=denominator > 0)


def holm_adjust(p_values):
    """Adjust the complete planned family, retaining undefined tests as p=1."""
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    adjusted = np.minimum(1, np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    result = np.empty_like(p)
    result[order] = adjusted
    return result


def paired_years(monthly, k, lead, month, flow_column="zero_fraction"):
    pairs = pd.DataFrame({"spi": monthly[f"SPI{k}"].shift(lead),
                          "flow": monthly[flow_column], "Q_n": monthly.Q_n})
    pairs = pairs.loc[START:END]
    pairs = pairs[(pairs.index.month == month) & pairs.Q_n.eq(pairs.index.days_in_month)].copy()
    pairs["HY"] = pairs.index.year - (pairs.index.month < 6)
    pairs = pairs.set_index("HY").reindex(YEARS)
    return pairs[["spi", "flow"]]


def correlation_analysis(monthly, rng, repetitions, permutations):
    samples = {length: block_indices(rng, 30, repetitions, length) for length in (1, 2, 3)}
    shuffle = np.array([rng.permutation(30) for _ in range(permutations)])
    block_order = np.array([rng.permutation(10) for _ in range(permutations)])
    block_shuffle = (block_order[:, :, None] * 3 + np.arange(3)).reshape(permutations, 30)
    rows = []
    for k in (3, 6):
        for lead in range(4):
            for month in range(1, 13):
                pair = paired_years(monthly, k, lead, month)
                x, y = pair.spi.to_numpy(), pair.flow.to_numpy()
                if not np.isfinite(pair.to_numpy()).all():
                    raise ValueError("Incomplete paired years need an explicit missing-data method.")
                row = dict(spi_scale=k, rainfall_lead_months=lead, flow_month=month,
                           n=30, distinct_flow_fractions=len(np.unique(y)), spearman_rho=np.nan,
                           permutation_p=np.nan, permutation_p_block3=np.nan,
                           constant_flow=bool(np.ptp(y) == 0))
                for length in (1, 2, 3):
                    row.update({f"ci_low_block{length}": np.nan, f"ci_high_block{length}": np.nan,
                                f"bootstrap_valid_fraction_block{length}": 0.0})
                if np.ptp(x) > 0 and np.ptp(y) > 0:
                    rho = float(stats.spearmanr(x, y).statistic)
                    row["spearman_rho"] = rho
                    xr, yr = stats.rankdata(x), stats.rankdata(y)
                    xr, yr = xr - xr.mean(), yr - yr.mean()
                    null = np.sum(xr[shuffle] * yr, axis=1) / np.sqrt(np.sum(xr*xr) * np.sum(yr*yr))
                    row["permutation_p"] = float((1 + np.sum(np.abs(null) >= abs(rho) - 1e-12)) /
                                                  (permutations + 1))
                    block_null = np.sum(xr[block_shuffle] * yr, axis=1) / np.sqrt(np.sum(xr*xr) * np.sum(yr*yr))
                    row["permutation_p_block3"] = float((1 + np.sum(np.abs(block_null) >= abs(rho) - 1e-12)) /
                                                         (permutations + 1))
                    for length, indices in samples.items():
                        boot = row_spearman(x[indices], y[indices])
                        valid = np.isfinite(boot)
                        fraction = float(valid.mean())
                        row[f"bootstrap_valid_fraction_block{length}"] = fraction
                        # Sparse dry-season variation can make many resamples constant.
                        if fraction >= 0.95:
                            low, high = np.quantile(boot[valid], [0.025, 0.975])
                            row[f"ci_low_block{length}"] = float(low)
                            row[f"ci_high_block{length}"] = float(high)
                    row["flow_lag1_rho"] = (float(stats.spearmanr(y[:-1], y[1:]).statistic)
                                             if np.ptp(y[:-1]) > 0 and np.ptp(y[1:]) > 0 else np.nan)
                    row["spi_lag1_rho"] = float(stats.spearmanr(x[:-1], x[1:]).statistic)
                rows.append(row)
    results = pd.DataFrame(rows)
    results["holm_p_96"] = holm_adjust(results.permutation_p.fillna(1))
    results["holm_reject_005"] = results.holm_p_96 < 0.05
    results["holm_p_96_block3"] = holm_adjust(results.permutation_p_block3.fillna(1))
    results["holm_reject_005_block3"] = results.holm_p_96_block3 < 0.05
    return results, samples[1]


def gamma_ks_statistic(values, shape, scale):
    values = np.sort(values)
    cdf = stats.gamma.cdf(values, shape, scale=scale)
    n = len(values)
    return float(max(np.max(np.arange(1, n+1)/n-cdf), np.max(cdf-np.arange(n)/n)))


def gamma_diagnostics(monthly, fits, rng, simulations, fit_bootstrap):
    diagnostics, uncertainty = [], []
    # The same sampled reference years preserve dependence between calendar-month fits.
    sampled_years = rng.integers(0, 40, size=(fit_bootstrap, 40))
    for fit in fits.itertuples(index=False):
        k, month = fit.spi_scale, fit.ending_month
        total = monthly[f"R{k}_mm"]
        ref = total[(total.index.year >= 1981) & (total.index.year <= 2020) &
                    (total.index.month == month)].to_numpy()
        pos = ref[ref > 0]
        observed_d = gamma_ks_statistic(pos, fit.shape, fit.scale_mm)
        simulated_d = []
        for _ in range(simulations):
            simulated = rng.gamma(fit.shape, fit.scale_mm, len(pos))
            shape, scale = gamma_parameters(simulated)
            simulated_d.append(gamma_ks_statistic(simulated, shape, scale))
        p_value = float((1 + np.sum(np.asarray(simulated_d) >= observed_d)) / (simulations + 1))
        diagnostics.append(dict(spi_scale=k, ending_month=month, n_reference=len(ref),
                                positive_n=len(pos), zero_probability=fit.zero_probability,
                                ks_positive_d=observed_d, refitted_mc_p=p_value))
        study = total.loc[START:END]
        study = study[study.index.month == month]
        transformed = np.full((fit_bootstrap, len(study)), np.nan)
        for i, sample in enumerate(ref[sampled_years]):
            if np.sum(sample > 0) >= 3 and np.ptp(sample[sample > 0]) > 0:
                shape, scale = gamma_parameters(sample)
                transformed[i] = gamma_spi(study.to_numpy(), shape, scale, np.mean(sample == 0))
        for j, date in enumerate(study.index):
            values = transformed[:, j]
            valid = np.isfinite(values)
            low, high = np.quantile(values[valid], [0.025, 0.975])
            uncertainty.append(dict(month=date.strftime("%Y-%m-%d"), spi_scale=k,
                                    spi=float(monthly.loc[date, f"SPI{k}"]),
                                    refit_ci_low=float(low), refit_ci_high=float(high),
                                    valid_fraction=float(valid.mean())))
    diagnostics = pd.DataFrame(diagnostics)
    diagnostics["holm_p_24"] = holm_adjust(diagnostics.refitted_mc_p)
    return diagnostics, pd.DataFrame(uncertainty).sort_values(["spi_scale", "month"])


def zero_spells(series):
    rows, first, last, length = [], None, None, 0
    for date, value in series.items():
        if pd.notna(value) and value == 0:
            first = date if first is None else first
            last, length = date, length + 1
        else:
            if length:
                rows.append(dict(start=first.date().isoformat(), end=last.date().isoformat(), days=length))
            first, last, length = None, None, 0
    if length:
        rows.append(dict(start=first.date().isoformat(), end=last.date().isoformat(), days=length))
    return pd.DataFrame(rows, columns=["start", "end", "days"]).sort_values(
        ["days", "start"], ascending=[False, True]).reset_index(drop=True)


def descriptive_tables(daily, bootstrap):
    study = daily.loc[START:END].copy()
    study["HY"] = study.index.year - (study.index.month < 6)
    annual = []
    seasonal = []
    climatology = []
    for year, group in study.groupby("HY"):
        runs = zero_spells(group.Q_m3s)
        annual.append(dict(HY=int(year), rainfall_mm=float(group.P_mm.sum()), days=len(group),
                           zero_days=int(group.Q_m3s.eq(0).sum()),
                           zero_percent=float(100 * group.Q_m3s.eq(0).mean()),
                           discharge_volume_Mm3=float(group.Q_m3s.sum() * 0.0864),
                           longest_zero_spell_within_HY_days=int(runs.days.max()) if len(runs) else 0))
    for name, months in SEASONS.items():
        group = study[study.index.month.isin(months)]
        by_year = group.groupby("HY").agg(days=("Q_m3s", "count"),
                                           zeros=("Q_m3s", lambda x: x.eq(0).sum()),
                                           volume=("Q_m3s", lambda x: x.sum() * 0.0864))
        bzero = 100 * by_year.zeros.to_numpy()[bootstrap].sum(axis=1) / by_year.days.to_numpy()[bootstrap].sum(axis=1)
        bvolume = by_year.volume.to_numpy()[bootstrap].mean(axis=1)
        seasonal.append(dict(season=name, mean_rainfall_mm=float(group.P_mm.sum()/30),
                             zero_percent=float(100*group.Q_m3s.eq(0).mean()),
                             zero_ci_low=float(np.quantile(bzero, .025)),
                             zero_ci_high=float(np.quantile(bzero, .975)),
                             mean_annual_season_volume_Mm3=float(by_year.volume.mean()),
                             volume_ci_low=float(np.quantile(bvolume, .025)),
                             volume_ci_high=float(np.quantile(bvolume, .975))))
    for month in MONTHS:
        group = study[study.index.month == month]
        by_year = group.groupby("HY").agg(days=("Q_m3s", "count"),
                                           zeros=("Q_m3s", lambda x: x.eq(0).sum()))
        boot = 100*by_year.zeros.to_numpy()[bootstrap].sum(axis=1)/by_year.days.to_numpy()[bootstrap].sum(axis=1)
        climatology.append(dict(month=month, rainfall_mm=float(group.P_mm.sum()/30),
                                zero_percent=float(100*group.Q_m3s.eq(0).mean()),
                                zero_ci_low=float(np.quantile(boot, .025)),
                                zero_ci_high=float(np.quantile(boot, .975))))
    return pd.DataFrame(annual), pd.DataFrame(seasonal), pd.DataFrame(climatology), zero_spells(study.Q_m3s)


def sensitivity_analysis(monthly, baseline):
    tables, events = [], []
    for name, reference, method in [("baseline", (1981, 2020), "mle"),
                                     ("reference_1981_2010", (1981, 2010), "mle"),
                                     ("gamma_moments", (1981, 2020), "moments")]:
        alternative, _ = calculate_spi(monthly, reference, method)
        for k in (3, 6):
            series = alternative[f"SPI{k}"].loc[START:END]
            episodes = drought_events(alternative[f"SPI{k}"])
            episodes.insert(0, "spi_scale", k)
            episodes.insert(0, "variant", name)
            events.append(episodes)
            tables.append(dict(variant=name, spi_scale=k, drought_episodes=len(episodes),
                               months_spi_le_minus1=int(series.le(-1).sum()),
                               minimum_spi=float(series.min()), peak_month=series.idxmin().strftime("%Y-%m"),
                               september_2003_spi=float(series.loc["2003-09-01"])))
    correlations = []
    for row in baseline.itertuples(index=False):
        pair = paired_years(monthly, row.spi_scale, row.rainfall_lead_months,
                            row.flow_month, "nearzero_fraction")
        rho = (float(stats.spearmanr(pair.spi, pair.flow).statistic)
               if pair.spi.nunique() > 1 and pair.flow.nunique() > 1 else np.nan)
        correlations.append(dict(spi_scale=row.spi_scale, rainfall_lead_months=row.rainfall_lead_months,
                                 flow_month=row.flow_month, zero_rho=row.spearman_rho,
                                 nearzero_rho=rho, difference=rho-row.spearman_rho))
    return pd.DataFrame(tables), pd.concat(events, ignore_index=True), pd.DataFrame(correlations)


def make_figures(monthly, climatology, correlations, spells, out):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "normal", "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(11.5, 7.5), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.15])
    rain, flow, timeline = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :])
    x = np.arange(12)
    rain.bar(x, climatology.rainfall_mm, color="#256880", width=.72)
    rain.set(title="(a) Mean monthly catchment rainfall", ylabel="Rainfall (mm/month)")
    flow.bar(x, climatology.zero_percent, color="#b46735", width=.72)
    flow.errorbar(x, climatology.zero_percent,
                  yerr=np.array([climatology.zero_percent-climatology.zero_ci_low,
                                 climatology.zero_ci_high-climatology.zero_percent]),
                  fmt="none", ecolor="#282828", capsize=3, linewidth=1)
    flow.set(title="(b) Recorded zero-flow frequency", ylabel="Zero-flow days (%)", ylim=(0, 109))
    for ax in (rain, flow):
        ax.set_xticks(x, LABELS, fontsize=10)
        ax.set_xlabel("Calendar month")
        ax.grid(axis="y", alpha=.2)
        ax.set_axisbelow(True)
    study = monthly.loc[START:END]
    timeline.axhline(0, color=".45", linewidth=.7)
    timeline.axhline(-1, color="#a84837", linestyle="--", linewidth=1, label="Drought intensity threshold (-1)")
    timeline.plot(study.index, study.SPI6, color="#256880", linewidth=1.1, label="SPI-6")
    timeline.fill_between(study.index, study.SPI6, 0, where=study.SPI6.lt(0), color="#b46735", alpha=.22)
    longest = spells.iloc[0]
    timeline.axvspan(pd.Timestamp(longest.start), pd.Timestamp(longest.end), color=".2", alpha=.15,
                     label=f"Longest zero-flow spell ({longest.days} days)")
    timeline.set(title="(c) Six-month rainfall deficits", ylabel="SPI-6 (dimensionless)", xlabel="Year",
                 xlim=(study.index[0], pd.Timestamp(END)), ylim=(-3.3, 3.3))
    timeline.xaxis.set_major_locator(mdates.YearLocator(5))
    timeline.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    timeline.legend(loc="upper left", ncol=3, frameon=False, fontsize=9)
    timeline.grid(axis="y", alpha=.2)
    fig.suptitle("Figure 1  Rainfall deficits and recorded flow intermittency at Sarati\nJune 1981 to May 2011", fontsize=15)
    fig.supxlabel("Whiskers in (b): 95% intervals from 10,000 resamples of complete hydrological years.  Source: CAMELS-IND v2.2, station 04049.", fontsize=9)
    fig.savefig(out / "Figure_1_Rainfall_and_Intermittency.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 8), layout="constrained", sharey=True)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e1e1e1")
    for ax, k in zip(axes, (3, 6)):
        sub = correlations[correlations.spi_scale.eq(k)]
        matrix = sub.pivot(index="flow_month", columns="rainfall_lead_months", values="spearman_rho").reindex(MONTHS)
        im = ax.imshow(matrix, vmin=-1, vmax=1, cmap=cmap, aspect="auto")
        ax.set(title=f"SPI-{k}", xlabel="Rainfall lead (months)")
        ax.set_xticks(range(4), [0, 1, 2, 3])
        ax.set_yticks(range(12), LABELS)
        ax.tick_params(labelleft=True)
        for i, month in enumerate(MONTHS):
            for lead in range(4):
                row = sub[sub.flow_month.eq(month) & sub.rainfall_lead_months.eq(lead)].iloc[0]
                text = "NA" if pd.isna(row.spearman_rho) else f"{row.spearman_rho:.2f}" + ("*" if row.holm_reject_005 else "")
                ax.text(lead, i, text, ha="center", va="center", fontsize=11,
                        color="white" if abs(row.spearman_rho) > .55 else "#202020")
        ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
        ax.set_yticks(np.arange(-.5, 12, 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1)
        ax.tick_params(which="minor", length=0)
    axes[0].set_ylabel("Month of discharge observation")
    fig.colorbar(im, ax=axes, label="Spearman correlation (rho)", shrink=.8, pad=.03)
    fig.suptitle("Figure 2  Rainfall SPI and monthly zero-flow fraction\n30 paired hydrological years per cell", fontsize=15)
    fig.supxlabel("* Holm-adjusted p < 0.05: independent-year permutation, 96 planned tests.  NA: constant flow fraction.\nSource: CAMELS-IND v2.2, station 04049. June 1981 to May 2011. See text for dependence sensitivity.", fontsize=8.5)
    fig.savefig(out / "Figure_2_Lagged_Associations.png", dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bootstrap", type=int, default=10000)
    parser.add_argument("--permutations", type=int, default=19999)
    parser.add_argument("--fit-simulations", type=int, default=999)
    parser.add_argument("--fit-bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=7662026)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    args = parser.parse_args()
    if min(args.bootstrap, args.permutations, args.fit_simulations, args.fit_bootstrap) < 99:
        parser.error("Use at least 99 repetitions for each calculation.")
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    seeds = np.random.SeedSequence(args.seed).spawn(2)
    daily = read_daily()
    monthly, fits = calculate_spi(monthly_data(daily))
    print("Calculating year-block intervals and permutation tests...", flush=True)
    correlations, bootstrap = correlation_analysis(monthly, np.random.default_rng(seeds[0]),
                                                   args.bootstrap, args.permutations)
    annual, seasonal, climatology, spells = descriptive_tables(daily, bootstrap)
    print("Checking gamma fits and reference-sample uncertainty...", flush=True)
    diagnostics, spi_uncertainty = gamma_diagnostics(monthly, fits, np.random.default_rng(seeds[1]),
                                                    args.fit_simulations, args.fit_bootstrap)
    sensitivity, events, threshold = sensitivity_analysis(monthly, correlations)
    study = daily.loc[START:END]
    baseline_events = events[events.variant.eq("baseline")]
    strongest = correlations.loc[correlations.spearman_rho.idxmin()].to_dict()
    summary = dict(analysis_first=START, analysis_last=END, hydrological_years=30,
                   n_days=len(study), zero_days=int(study.Q_m3s.eq(0).sum()),
                   zero_percent=float(100*study.Q_m3s.eq(0).mean()),
                   nearzero_percent=float(100*study.Q_m3s.le(.1).mean()),
                   mean_annual_rainfall_mm=float(study.P_mm.sum()/30),
                   monsoon_rainfall_share_percent=float(100*study.loc[study.index.month.isin(SEASONS["Jun-Sep"]), "P_mm"].sum()/study.P_mm.sum()),
                   longest_zero_spell=spells.iloc[0].to_dict(),
                   strongest_negative_association=strongest,
                   planned_tests=96, defined_tests=int(correlations.spearman_rho.notna().sum()),
                   holm_significant_tests=int(correlations.holm_reject_005.sum()),
                   holm_significant_tests_block3=int(correlations.holm_reject_005_block3.sum()),
                   gamma_fits_rejected_unadjusted=int(diagnostics.refitted_mc_p.lt(.05).sum()),
                   gamma_fits_rejected_holm=int(diagnostics.holm_p_24.lt(.05).sum()),
                   nearzero_max_absolute_rho_change=float(threshold.difference.abs().max()),
                   spi3_events=int(baseline_events.spi_scale.eq(3).sum()),
                   spi6_events=int(baseline_events.spi_scale.eq(6).sum()),
                   spi6_largest_magnitude_event=baseline_events[baseline_events.spi_scale.eq(6)].sort_values("magnitude", ascending=False).iloc[0].to_dict(),
                   min_positive_discharge_m3s=float(study.loc[study.Q_m3s.gt(0), "Q_m3s"].min()),
                   input_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in sorted(RAW.rglob("*.csv"))})
    config = dict(seed=args.seed, bootstrap_repetitions=args.bootstrap,
                  year_block_lengths=[1, 2, 3], permutation_repetitions=args.permutations,
                  gamma_fit_mc_repetitions=args.fit_simulations,
                  spi_reference_bootstrap_repetitions=args.fit_bootstrap,
                  python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                  scipy=scipy.__version__, matplotlib=matplotlib.__version__)
    exports = {"monthly_spi_and_flow": monthly, "gamma_parameters": fits,
               "lagged_correlations": correlations, "annual_metrics": annual,
               "seasonal_metrics": seasonal, "monthly_climatology": climatology,
               "zero_flow_spells": spells, "gamma_fit_diagnostics": diagnostics,
               "spi_reference_uncertainty": spi_uncertainty, "spi_sensitivity": sensitivity,
               "drought_events": events, "nearzero_correlation_sensitivity": threshold}
    for name, frame in exports.items():
        frame.to_csv(out / f"{name}.csv", index=(name == "monthly_spi_and_flow"),
                     index_label="month" if name == "monthly_spi_and_flow" else None)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    (out / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    make_figures(monthly, climatology, correlations, spells, out)
    print(json.dumps({key: value for key, value in summary.items() if key != "input_sha256"}, indent=2))


if __name__ == "__main__":
    main()
