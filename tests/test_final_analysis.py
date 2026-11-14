"""Checks for tied ranks, year boundaries and drought-event definitions."""
# Author: Ganesh Sonawane
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
import final_analysis as final


class StatisticalChecks(unittest.TestCase):
    def test_gamma_ks_statistic_matches_scipy(self):
        values = np.array([.02, .02, .4, 1, 4.5, 9])
        expected = stats.kstest(values, "gamma", args=(.8, 0, 3)).statistic
        self.assertAlmostEqual(final.gamma_ks_statistic(values, .8, 3), expected, places=14)

    def test_bootstrap_repeated_values_are_ranked_again(self):
        x = np.array([[1, 1, 4, 8, 9], [9, 4, 4, 1, 1]], dtype=float)
        y = np.array([[0, 0, 1, 0, 2], [2, 1, 1, 0, 0]], dtype=float)
        expected = [stats.spearmanr(a, b).statistic for a, b in zip(x, y)]
        np.testing.assert_allclose(final.row_spearman(x, y), expected, atol=1e-14)

    def test_constant_resamples_stay_undefined(self):
        x = np.array([[1, 2, 3], [1, 1, 1]], dtype=float)
        y = np.array([[0, 0, 0], [1, 2, 3]], dtype=float)
        self.assertTrue(np.isnan(final.row_spearman(x, y)).all())

    def test_holm_keeps_family_size_and_monotonicity(self):
        np.testing.assert_allclose(final.holm_adjust([.5, .001, .02]), [.5, .003, .04])
        p = np.ones(96)
        p[7] = .0001
        self.assertAlmostEqual(final.holm_adjust(p)[7], .0096)

    def test_circular_blocks_keep_adjacent_years(self):
        samples = final.block_indices(np.random.default_rng(71), 30, 100, 3)
        self.assertEqual(samples.shape, (100, 30))
        blocks = samples.reshape(100, 10, 3)
        self.assertTrue(np.all(np.diff(blocks, axis=2) % 30 == 1))
        self.assertTrue(np.all((samples >= 0) & (samples < 30)))


class EventChecks(unittest.TestCase):
    def test_negative_lead_in_counts_but_weak_spell_does_not(self):
        values = pd.Series([.3, -.2, -1.4, -.1, .2, -.6, -.4, .1],
                           index=pd.date_range("2000-01-01", periods=8, freq="MS"))
        events = final.drought_events(values, "2000-01-01", "2000-08-01")
        self.assertEqual(len(events), 1)
        row = events.iloc[0]
        self.assertEqual((row.start, row.end, row.duration_months), ("2000-02", "2000-04", 3))
        self.assertAlmostEqual(row.magnitude, 1.7)
        self.assertEqual(row.peak_month, "2000-03")
        self.assertFalse(row.left_censored or row.right_censored)

    def test_boundary_events_are_flagged(self):
        values = pd.Series([-.2, -1.2, -.4, .2, -.2, -1.4, -.7],
                           index=pd.date_range("2000-01-01", periods=7, freq="MS"))
        events = final.drought_events(values, "2000-02-01", "2000-06-01")
        self.assertEqual(len(events), 2)
        self.assertTrue(events.iloc[0].left_censored)
        self.assertTrue(events.iloc[1].right_censored)

    def test_zero_or_missing_spi_ends_a_negative_sequence(self):
        values = pd.Series([-1.2, 0, -.4, np.nan, -.7, -1.1, .4],
                           index=pd.date_range("2000-01-01", periods=7, freq="MS"))
        events = final.drought_events(values, "2000-01-01", "2000-07-01")
        self.assertEqual(list(events.duration_months), [1, 2])

    def test_recorded_zero_spells_include_leap_day_and_break_on_missing(self):
        values = pd.Series([0, 0, 0, np.nan, 0, .01],
                           index=pd.date_range("2004-02-28", periods=6))
        spells = final.zero_spells(values)
        self.assertEqual(list(spells.days), [3, 1])
        self.assertEqual(spells.iloc[0].end, "2004-03-01")
        self.assertTrue(final.zero_spells(values + 1).empty)


class StudyChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.monthly, _ = final.calculate_spi(final.monthly_data(final.read_daily()))

    def test_lags_are_aligned_before_year_resampling(self):
        pair = final.paired_years(self.monthly, 6, 3, 6)
        self.assertEqual(len(pair), 30)
        self.assertAlmostEqual(pair.loc[1981, "spi"], self.monthly.loc["1981-03-01", "SPI6"])
        self.assertAlmostEqual(pair.loc[2010, "flow"], self.monthly.loc["2010-06-01", "zero_fraction"])
        january = final.paired_years(self.monthly, 3, 1, 1)
        self.assertAlmostEqual(january.loc[1981, "spi"], self.monthly.loc["1981-12-01", "SPI3"])

    def test_preliminary_spi_values_are_preserved(self):
        previous = pd.read_csv(final.ROOT / "analysis/results/monthly_spi_and_flow.csv",
                               index_col=0, parse_dates=True)
        for column in ("SPI3", "SPI6", "zero_fraction"):
            np.testing.assert_allclose(self.monthly[column], previous[column], atol=1e-12, equal_nan=True)

    def test_reference_change_keeps_month_specific_rank_order(self):
        shorter, _ = final.calculate_spi(self.monthly, (1981, 2010))
        for k in (3, 6):
            for month in range(1, 13):
                original = final.paired_years(self.monthly, k, 2, month)
                alternative = final.paired_years(shorter, k, 2, month)
                np.testing.assert_array_equal(stats.rankdata(original.spi), stats.rankdata(alternative.spi))

    def test_identical_dry_window_totals_keep_identical_spi(self):
        february = self.monthly[self.monthly.index.month == 2]
        tied = february[february.R3_mm.eq(.02)]
        self.assertGreater(len(tied), 1)
        self.assertEqual(tied.SPI3.nunique(), 1)


if __name__ == "__main__":
    unittest.main()
