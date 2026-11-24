"""Route measured Sarati flow through a hypothetical seasonal storage scheme."""
# Author: Ganesh Sonawane
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from itertools import product
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis/water_budget_results"
DAY = 86400.0
AREA_KM2 = 6807.1
START, END = "1981-06-01", "2011-05-31"
FLOW_NAMES = ["inflow", "diversion", "rain", "evaporation", "seepage", "release",
              "transfer_loss", "delivered", "spill", "bypass", "downstream"]


@dataclass(frozen=True)
class Design:
    capacity_Mm3: float = 10.0
    capture_fraction: float = 0.10
    target_m3s: float = 1.0
    bypass_floor_m3s: float = 1.0
    intake_limit_m3s: float = 5.0
    release_limit_m3s: float = 2.0
    full_pool_mean_depth_m: float = 10.0
    area_exponent: float = 2/3
    seepage_fraction_per_day: float = 0.001
    delivery_efficiency: float = 0.9
    evaporation_factor: float = 1.0
    pet_source: str = "Singer"
    initial_fraction: float = 0.0


def read_inputs():
    rain = pd.read_csv(ROOT / "data/raw/water_budget_forcings/04049.csv")
    flow = pd.read_csv(ROOT / "data/raw/sarati_source/streamflow_timeseries/streamflow_observed.csv")
    for frame in (rain, flow):
        frame.index = pd.to_datetime(frame[["year", "month", "day"]])
        if not frame.index.equals(pd.date_range("1980-01-01", "2020-12-31")):
            raise ValueError("Input dates are duplicated, missing or out of order.")
    data = pd.DataFrame({"rain_mm": rain["prcp(mm/day)"], "pet_Singer_mm": rain["pet(mm/day)"],
                         "pet_GLEAM_mm": rain["pet_gleam(mm/day)"],
                         "aet_GLEAM_mm": rain["aet_gleam(mm/day)"], "Q_m3s": flow["4049"]}).loc[START:END]
    if len(data) != 10957 or not np.isfinite(data.to_numpy()).all() or (data < 0).any().any():
        raise ValueError("The model needs 10,957 nonnegative, complete input days.")
    return data


def step(storage, q, rain_mm, pet_mm, monsoon, design, dt=1.0):
    """Volumes are m3. Daily means are held constant within each substep."""
    d = design
    capacity = d.capacity_Mm3 * 1e6
    if not (0 <= storage <= capacity + 1e-6 and 0 < dt <= 1):
        raise ValueError("Invalid storage or timestep.")
    inflow = q * DAY * dt
    if capacity == 0:
        return 0.0, dict(zip(FLOW_NAMES, [inflow, 0, 0, 0, 0, 0, 0, 0, 0, inflow, inflow]))
    potential = d.capture_fraction * max(q - d.bypass_floor_m3s, 0) * DAY * dt
    diversion = min(potential, d.intake_limit_m3s * DAY * dt, capacity-storage) if monsoon else 0.0
    wet_storage = storage + diversion
    area = capacity/d.full_pool_mean_depth_m * (wet_storage/capacity)**d.area_exponent
    rain = rain_mm/1000 * area * dt
    evaporation = min(d.evaporation_factor * pet_mm/1000 * area * dt, wet_storage+rain)
    available = wet_storage + rain - evaporation
    seepage = min(available, available * (1-(1-d.seepage_fraction_per_day)**dt))
    available -= seepage
    requested = max(d.target_m3s-q, 0) * DAY * dt / d.delivery_efficiency
    release = min(requested, d.release_limit_m3s * DAY * dt, available) if not monsoon else 0.0
    available -= release
    spill = max(available-capacity, 0)
    end = available-spill
    delivered = d.delivery_efficiency * release
    transfer_loss = release-delivered
    bypass = inflow-diversion
    downstream = bypass + delivered + spill
    values = [inflow, diversion, rain, evaporation, seepage, release, transfer_loss,
              delivered, spill, bypass, downstream]
    return end, dict(zip(FLOW_NAMES, values))


def simulate(data, design=Design(), dt=1.0):
    d = design
    if (d.capacity_Mm3 < 0 or not 0 <= d.capture_fraction <= 1 or
            not 0 < d.delivery_efficiency <= 1 or not 0 <= d.initial_fraction <= 1 or
            not 0 <= d.seepage_fraction_per_day < 1 or d.full_pool_mean_depth_m <= 0):
        raise ValueError("Invalid design parameters.")
    substeps = round(1/dt)
    if abs(dt*substeps-1) > 1e-12:
        raise ValueError("Substeps must divide one day exactly.")
    storage = d.capacity_Mm3 * 1e6 * d.initial_fraction
    output = []
    fields = data[["Q_m3s", "rain_mm", f"pet_{d.pet_source}_mm"]].to_numpy()
    for date, (q, rain, pet) in zip(data.index, fields):
        initial = storage
        totals = {name: 0.0 for name in FLOW_NAMES}
        for _ in range(substeps):
            storage, flux = step(storage, q, rain, pet, date.month in (6, 7, 8, 9), d, dt)
            for name, value in flux.items():
                totals[name] += value
        reservoir_residual = (storage-initial)-(totals["diversion"]+totals["rain"]-
                                                totals["evaporation"]-totals["seepage"]-
                                                totals["release"]-totals["spill"])
        system_residual = (totals["inflow"]+totals["rain"]-totals["downstream"]-
                           totals["evaporation"]-totals["seepage"]-
                           totals["transfer_loss"]-(storage-initial))
        output.append([initial, storage, *totals.values(), reservoir_residual, system_residual])
    columns = ["storage_start", "storage_end", *FLOW_NAMES, "reservoir_residual", "system_residual"]
    result = pd.DataFrame(output, index=data.index, columns=[name+"_m3" for name in columns])
    result["Q_downstream_m3s"] = result.downstream_m3/DAY
    if (result.storage_end_m3.min() < -1e-7 or result.storage_end_m3.max() > d.capacity_Mm3*1e6+1e-7 or
            result.reservoir_residual_m3.abs().max() > 1e-5 or result.system_residual_m3.abs().max() > 1e-5):
        raise ArithmeticError("Mass balance or storage constraint failed.")
    return result


def metrics(data, run, design):
    dry = ~data.index.month.isin([6, 7, 8, 9])
    monsoon = ~dry
    target = design.target_m3s
    deficit = np.maximum(target-run.Q_downstream_m3s.to_numpy(), 0)*DAY
    base_deficit = np.maximum(target-data.Q_m3s.to_numpy(), 0)*DAY
    sums = {name+"_Mm3_per_year": float(run[name+"_m3"].sum()/30e6) for name in FLOW_NAMES}
    sums.update(dry_deficit_Mm3_per_year=float(deficit[dry].sum()/30e6),
                baseline_dry_deficit_Mm3_per_year=float(base_deficit[dry].sum()/30e6),
                dry_deficit_reduction_percent=float(100*(1-deficit[dry].sum()/base_deficit[dry].sum())),
                dry_target_reliability_percent=float(100*(run.Q_downstream_m3s[dry] >= target-1e-9).mean()),
                dry_zero_percent=float(100*(run.Q_downstream_m3s[dry] <= 1e-10).mean()),
                monsoon_diversion_percent=float(100*run.diversion_m3[monsoon].sum()/run.inflow_m3[monsoon].sum()),
                terminal_storage_Mm3=float(run.storage_end_m3.iloc[-1]/1e6),
                initial_storage_Mm3=float(run.storage_start_m3.iloc[0]/1e6),
                storage_change_Mm3_per_year=float((run.storage_end_m3.iloc[-1]-run.storage_start_m3.iloc[0])/30e6),
                max_reservoir_residual_m3=float(run.reservoir_residual_m3.abs().max()),
                max_system_residual_m3=float(run.system_residual_m3.abs().max()))
    return sums


def annual_ledger(data, run, name):
    groups = data.index.year-(data.index.month < 6)
    rows = []
    for hy in np.unique(groups):
        sub = run[groups == hy]
        dry = ~sub.index.month.isin([6, 7, 8, 9])
        row = {"scenario": name, "HY": int(hy)}
        row.update({term+"_Mm3": float(sub[term+"_m3"].sum()/1e6) for term in FLOW_NAMES})
        row.update(storage_start_Mm3=float(sub.storage_start_m3.iloc[0]/1e6),
                   storage_end_Mm3=float(sub.storage_end_m3.iloc[-1]/1e6),
                   dry_deficit_Mm3=float(np.maximum(1-sub.Q_downstream_m3s[dry], 0).sum()*DAY/1e6),
                   reservoir_residual_m3=float(sub.reservoir_residual_m3.sum()),
                   system_residual_m3=float(sub.system_residual_m3.sum()))
        rows.append(row)
    return pd.DataFrame(rows)


def watershed_diagnostic(data):
    hy = data.index.year-(data.index.month < 6)
    rows = []
    for year in np.unique(hy):
        sub = data[hy == year]
        p, et = sub.rain_mm.sum(), sub.aet_GLEAM_mm.sum()
        runoff = sub.Q_m3s.sum()*DAY/(AREA_KM2*1000)
        rows.append(dict(HY=int(year), rainfall_mm=p, aet_GLEAM_mm=et, observed_runoff_mm=runoff,
                         unclosed_residual_mm=p-et-runoff))
    return pd.DataFrame(rows)


def figures(data, nominal, ranges, runs, out):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.1), layout="constrained")
    terms = [("delivered", "Delivered return", "#2d7893"), ("transfer_loss", "Transfer loss", "#80b1c4"),
             ("evaporation", "Evaporation", "#cb8b45"), ("seepage", "Seepage", "#c85d45"),
             ("spill", "Spill", "#999999"), ("storage_change", "Net storage change", "#594d7a")]
    middle = nominal[nominal.capture_fraction.eq(.1)].sort_values("capacity_Mm3")
    centers = np.arange(3)
    incoming = middle.diversion_Mm3_per_year+middle.rain_Mm3_per_year
    axes[0].bar(centers-.19, incoming, width=.36, color="#326647", label="Diversion + direct rain")
    bottom = np.zeros(3)
    for key, label, color in terms:
        value = middle[key+"_Mm3_per_year"].to_numpy()
        axes[0].bar(centers+.19, value, bottom=bottom, width=.36, color=color, label=label)
        bottom += value
    axes[0].set_xticks(centers, ["5", "10", "20"])
    axes[0].set(xlabel="Assumed storage capacity (million m3)", ylabel="Mean volume (million m3/year)",
                title="(a) Storage budget at 10% capture")
    axes[0].legend(fontsize=8.5, frameon=False, loc="upper left")
    for i, row in enumerate(middle.itertuples(index=False)):
        axes[0].text(i-.19, incoming.iloc[i]+.35, f"{incoming.iloc[i]:.2f}", ha="center", fontsize=9)
    colors = ["#7d6a8c", "#277f96", "#c17839"]
    for color, fraction in zip(colors, [.05, .1, .2]):
        sub = nominal[nominal.capture_fraction.eq(fraction)].sort_values("capacity_Mm3")
        band = ranges[ranges.capture_fraction.eq(fraction)].sort_values("capacity_Mm3")
        x = sub.capacity_Mm3.to_numpy()
        axes[1].plot(x, sub.dry_deficit_reduction_percent, "o-", color=color, label=f"Capture {100*fraction:.0f}%")
        axes[1].fill_between(x, band.reduction_min, band.reduction_max, color=color, alpha=.13)
    axes[1].set(xlabel="Assumed storage capacity (million m3)", ylabel="Oct-May deficit-volume reduction (%)",
                title="(b) Intervention performance", ylim=(0, 105), xticks=[5, 10, 20])
    axes[1].legend(frameon=False, fontsize=9)
    for ax in axes:
        ax.grid(axis="y", alpha=.18)
        ax.set_axisbelow(True)
    fig.suptitle("Figure 1  Water budget and seasonal storage interventions at Sarati", fontsize=15)
    fig.supxlabel("Hypothetical scheme, June 1981-May 2011. Target 1 m3/s. Bands: 36 assumed loss/PET combinations, not confidence intervals.", fontsize=9)
    fig.savefig(out/"Figure_1_Water_Budget_and_Interventions.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(11.5, 8), layout="constrained", sharex=True)
    a, b = "2002-05-01", "2004-04-30"
    subset = data.loc[a:b]
    for color, capacity in zip(colors, [5, 10, 20]):
        run = runs[f"K{capacity}_f10"].loc[a:b]
        axes[0].plot(run.index, run.Q_downstream_m3s, color=color, linewidth=.9, alpha=.8, label=f"{capacity} million m3")
        axes[1].plot(run.index, run.storage_end_m3/1e6, color=color, label=f"{capacity} million m3")
    axes[0].plot(subset.index, subset.Q_m3s, color=".15", linestyle="--", linewidth=1.3,
                 label="Measured baseline", zorder=5)
    axes[0].axhline(1, color=".4", linestyle="--", linewidth=.8)
    axes[0].set(ylabel="Flow (m3/s)", title="(a) Baseline and simulated downstream flow", ylim=(0, 2.2))
    axes[0].legend(ncol=4, frameon=False, fontsize=9)
    axes[1].set(ylabel="Storage (million m3)", title="(b) Simulated storage at 10% capture")
    axes[1].text(.99, .90, "Capacity curves nearly overlap: inflow limits refill.\nAll three stores empty on 3 Nov 2002 and 7 Oct 2003.",
                 transform=axes[1].transAxes, ha="right", va="top", fontsize=9)
    axes[2].plot(subset.index, subset.rain_mm, color="#326647", linewidth=.8)
    axes[2].set(ylabel="Rainfall (mm/day)", xlabel="Date", title="(c) Catchment-mean rainfall")
    for ax in axes:
        for year in [2002, 2003]:
            ax.axvspan(pd.Timestamp(year, 6, 1), pd.Timestamp(year, 9, 30), color="#a7c9c9", alpha=.22)
        ax.grid(axis="y", alpha=.18)
        ax.set_xlim(pd.Timestamp(a), pd.Timestamp(b))
    fig.suptitle("Figure 2  Intervention behavior during the 2002-2004 rainfall drought", fontsize=15)
    fig.supxlabel("Central assumptions; empty start in June 1981 with continuous carryover. Shading: diversion season.\nFlow axis limited to 2.2 m3/s; higher flows are clipped. Source inputs: CAMELS-IND v2.2, station 04049.", fontsize=9)
    fig.savefig(out/"Figure_2_Drought_Intervention_Trajectories.png", dpi=300)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    data = read_inputs()
    baseline_design = Design(capacity_Mm3=0)
    baseline = simulate(data, baseline_design)
    baseline.to_csv(OUT/"daily_baseline.csv", index_label="date")
    baseline_metrics = metrics(data, baseline, baseline_design)
    nominal_rows, ensemble_rows, checks, ledgers = [], [], [], []
    runs = {}
    for capacity, fraction in product([5, 10, 20], [.05, .1, .2]):
        name = f"K{capacity}_f{round(fraction*100)}"
        d = Design(capacity_Mm3=capacity, capture_fraction=fraction)
        run = simulate(data, d)
        runs[name] = run
        run.to_csv(OUT/f"daily_{name}.csv", index_label="date")
        nominal_rows.append(dict(scenario=name, **asdict(d), **metrics(data, run, d)))
        ledgers.append(annual_ledger(data, run, name))
        for source_factor, seepage, eta in product([("Singer", .8), ("Singer", 1), ("Singer", 1.2), ("GLEAM", 1)],
                                                  [.0005, .001, .002], [.8, .9, 1]):
            source, factor = source_factor
            alternative = replace(d, pet_source=source, evaporation_factor=factor,
                                  seepage_fraction_per_day=seepage, delivery_efficiency=eta)
            result = run if alternative == d else simulate(data, alternative)
            ensemble_rows.append(dict(scenario=name, **asdict(alternative), **metrics(data, result, alternative)))
        for label, variant, dt in [("initial_full", replace(d, initial_fraction=1), 1),
                                   ("target_0.5", replace(d, target_m3s=.5), 1),
                                   ("target_2", replace(d, target_m3s=2), 1),
                                   ("half_day", d, .5)]:
            result = simulate(data, variant, dt)
            checks.append(dict(scenario=name, check=label, **asdict(variant), **metrics(data, result, variant)))
        print("Completed", name, flush=True)
    nominal = pd.DataFrame(nominal_rows)
    ensemble = pd.DataFrame(ensemble_rows)
    ranges = ensemble.groupby(["capacity_Mm3", "capture_fraction"], as_index=False).agg(
        reduction_min=("dry_deficit_reduction_percent", "min"), reduction_max=("dry_deficit_reduction_percent", "max"),
        reliability_min=("dry_target_reliability_percent", "min"), reliability_max=("dry_target_reliability_percent", "max"))
    sensitivity = pd.DataFrame(checks)
    diagnostic = watershed_diagnostic(data)
    nominal.to_csv(OUT/"nominal_designs.csv", index=False)
    ensemble.to_csv(OUT/"uncertainty_ensemble.csv", index=False)
    ranges.to_csv(OUT/"ensemble_ranges.csv", index=False)
    sensitivity.to_csv(OUT/"target_initial_timestep_checks.csv", index=False)
    pd.concat(ledgers, ignore_index=True).to_csv(OUT/"annual_storage_ledgers.csv", index=False)
    diagnostic.to_csv(OUT/"catchment_budget_diagnostic.csv", index=False)
    central = nominal[nominal.scenario.eq("K10_f10")].iloc[0].to_dict()
    source_paths = [ROOT/"data/raw/water_budget_forcings/04049.csv",
                    ROOT/"data/raw/sarati_source/streamflow_timeseries/streamflow_observed.csv"]
    summary = dict(period=[START, END], days=len(data), baseline=baseline_metrics, central_design=central,
                   ensemble_runs=len(ensemble), targeted_check_runs=len(sensitivity),
                   catchment_mean_budget_mm=diagnostic.drop(columns="HY").mean().to_dict(),
                   halfday_max_reduction_change_pp=float(np.max(np.abs(sensitivity[sensitivity.check.eq("half_day")].dry_deficit_reduction_percent.to_numpy()-nominal.dry_deficit_reduction_percent.to_numpy()))),
                   initial_full_max_reduction_change_pp=float(np.max(np.abs(sensitivity[sensitivity.check.eq("initial_full")].dry_deficit_reduction_percent.to_numpy()-nominal.dry_deficit_reduction_percent.to_numpy()))),
                   input_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths})
    (OUT/"summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    config = dict(central_assumptions=asdict(Design()), capture_season="June-September", release_season="October-May",
                  timestep_days=1, ensemble="36 equally unweighted combinations per design; min-max ranges have no probability interpretation",
                  control_volume="Diversion junction and new off-channel storage just downstream of the Sarati gauge",
                  source_flow="Observed regulated Sarati flow held fixed as upstream boundary; no extra catchment runoff is added",
                  step_order="diversion; wetted area; rain; evaporation; seepage; release; spill",
                  area_rule="A=(K/10m)*(S_after_diversion/K)^(2/3); evaluated once per step",
                  caveat="Uncalibrated hypothetical geometry, thresholds and losses; results are conditional screening, not a reservoir design")
    (OUT/"model_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    figures(data, nominal, ranges, runs, OUT)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
