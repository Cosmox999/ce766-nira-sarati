"""Cross-check final tables against raw CSV values and separate calculations."""
# Author: Ganesh Sonawane
import csv
import datetime as dt
import hashlib
import json
import math
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from statistics import NormalDist

import numpy as np
from scipy import stats
from scipy.special import digamma, gammainc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis/final_results"


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def month_offset(year, month, offset):
    serial = year * 12 + month - 1 + offset
    return serial // 12, serial % 12 + 1


def main():
    summary = json.loads((OUT / "summary.json").read_text())
    for path, checksum in summary["input_sha256"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == checksum, path
    rainfall = defaultdict(Decimal)
    flow = defaultdict(list)
    for row in read_csv(ROOT / "data/raw/sarati_source/catchment_mean_forcings/04049.csv"):
        rainfall[int(row["year"]), int(row["month"])] += Decimal(row["prcp(mm/day)"])
    for row in read_csv(ROOT / "data/raw/sarati_source/streamflow_timeseries/streamflow_observed.csv"):
        date = dt.date(int(row["year"]), int(row["month"]), int(row["day"]))
        if dt.date(1981, 6, 1) <= date <= dt.date(2011, 5, 31):
            flow[date.year, date.month].append(Decimal(row["4049"]))
    assert sum(map(len, flow.values())) == summary["n_days"] == 10957
    assert sum(value == 0 for values in flow.values() for value in values) == summary["zero_days"] == 6984
    assert abs(sum(value <= Decimal('.1') for values in flow.values() for value in values) * 100 / 10957
               - summary["nearzero_percent"]) < 1e-12
    monthly = {(int(row["month"][:4]), int(row["month"][5:7])): row
               for row in read_csv(OUT / "monthly_spi_and_flow.csv")}
    for fit in read_csv(OUT / "gamma_parameters.csv"):
        k, month = int(fit["spi_scale"]), int(fit["ending_month"])
        ref = [float(sum(rainfall[month_offset(year, month, -i)] for i in range(k)))
               for year in range(1981, 2021)]
        pos = [value for value in ref if value > 0]
        shape, scale, p0 = float(fit["shape"]), float(fit["scale_mm"]), float(fit["zero_probability"])
        assert abs(scale - np.mean(pos)/shape) < 1e-9
        score = math.log(shape) - digamma(shape) - math.log(np.mean(pos)) + np.mean(np.log(pos))
        assert abs(score) < 1e-8
        assert p0 == sum(value == 0 for value in ref)/40
        for (year, ending), values in flow.items():
            if ending != month:
                continue
            amount = float(sum(rainfall[month_offset(year, month, -i)] for i in range(k)))
            h = max(1e-10, min(1 - 1e-10, p0 + (1-p0)*gammainc(shape, amount/scale)))
            expected = NormalDist().inv_cdf(h)
            assert abs(expected - float(monthly[year, month][f"SPI{k}"])) < 1e-9
    for row in read_csv(OUT / "lagged_correlations.csv"):
        k, lead, month = int(row["spi_scale"]), int(row["rainfall_lead_months"]), int(row["flow_month"])
        x, y = [], []
        for hy in range(1981, 2011):
            year = hy + int(month < 6)
            x.append(float(sum(rainfall[month_offset(year, month, -lead-i)] for i in range(k))))
            values = flow[year, month]
            y.append(sum(value == 0 for value in values)/len(values))
        if len(set(y)) == 1:
            assert row["spearman_rho"] == ""
        else:
            expected = stats.spearmanr(x, y).statistic
            assert abs(expected - float(row["spearman_rho"])) < 1e-12
        for length in (1, 2, 3):
            valid = float(row[f"bootstrap_valid_fraction_block{length}"])
            if valid < .95:
                assert row[f"ci_low_block{length}"] == row[f"ci_high_block{length}"] == ""
            else:
                assert -1 <= float(row[f"ci_low_block{length}"]) <= float(row[f"ci_high_block{length}"]) <= 1
    for row in read_csv(OUT / "seasonal_metrics.csv"):
        months = {"Jun-Sep": [6, 7, 8, 9], "Oct-Nov": [10, 11],
                  "Dec-Feb": [12, 1, 2], "Mar-May": [3, 4, 5]}[row["season"]]
        values = [value for (year, month), days in flow.items() if month in months for value in days]
        volume = float(sum(values) * Decimal('0.0864')/30)
        assert abs(volume - float(row["mean_annual_season_volume_Mm3"])) < 1e-8
    events = read_csv(OUT / "drought_events.csv")
    for row in events:
        if row["variant"] != "baseline":
            continue
        start, end = tuple(map(int, row["start"].split("-"))), tuple(map(int, row["end"].split("-")))
        values = [float(record[f"SPI{row['spi_scale']}"]) for key, record in monthly.items() if start <= key <= end]
        assert len(values) == int(row["duration_months"])
        assert max(values) < 0 and min(values) <= -1
        assert abs(-sum(values) - float(row["magnitude"])) < 1e-10
    result = ("PASS: source checksums, raw counts, all study SPI values, all 96 rainfall-flow rank pairs, "
              "constant-series handling, interval availability, seasonal mean volumes and baseline drought episodes.\n")
    (OUT / "verification.txt").write_text(result, encoding="utf-8")
    print(result.strip())


if __name__ == "__main__":
    main()
