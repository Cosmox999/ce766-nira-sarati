"""Conservation and operating-rule checks for the storage model."""
# Author: Ganesh Sonawane
from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"analysis"))
from water_budget import Design, DAY, FLOW_NAMES, simulate, step


class BudgetChecks(unittest.TestCase):
    def assert_balanced(self, initial, q, final, flux):
        reservoir = flux['diversion']+flux['rain']-flux['evaporation']-flux['seepage']-flux['release']-flux['spill']
        self.assertAlmostEqual(final-initial, reservoir, places=7)
        system = flux['inflow']+flux['rain']-flux['downstream']-flux['evaporation']-flux['seepage']-flux['transfer_loss']
        self.assertAlmostEqual(final-initial, system, places=7)
        self.assertTrue(all(value >= 0 for value in flux.values()))

    def test_baseline_is_exact_pass_through(self):
        end, f = step(0, 27, 40, 8, True, Design(capacity_Mm3=0))
        self.assertEqual(end, 0)
        self.assertEqual(f['inflow'], f['downstream'])
        self.assertEqual(f['diversion']+f['rain']+f['release'], 0)

    def test_dry_empty_storage_cannot_create_release(self):
        end, f = step(0, 0, 0, 5, False, Design())
        self.assertEqual(end+f['release']+f['downstream'], 0)

    def test_diversion_never_uses_reserved_bypass(self):
        d = Design()
        end, f = step(0, .9, 0, 0, True, d)
        self.assertEqual(f['diversion'], 0)
        end, f = step(0, 11, 0, 0, True, replace(d, seepage_fraction_per_day=0))
        self.assertEqual(f['diversion'], DAY)
        self.assertGreaterEqual(f['bypass']/DAY, d.bypass_floor_m3s)
        self.assert_balanced(0, 11, end, f)

    def test_storage_and_intake_limits(self):
        d = Design(seepage_fraction_per_day=0)
        end, f = step(0, 1000, 0, 0, True, d)
        self.assertEqual(f['diversion'], d.intake_limit_m3s*DAY)
        end, f = step(10e6-100, 1000, 0, 0, True, d)
        self.assertEqual(f['diversion'], 100)
        self.assertEqual(end, 10e6)

    def test_lossless_release_meets_target_and_consumes_storage(self):
        d = Design(seepage_fraction_per_day=0, delivery_efficiency=1)
        end, f = step(1e6, 0, 0, 0, False, d)
        self.assertEqual(f['release'], DAY)
        self.assertEqual(end, 1e6-DAY)
        self.assertEqual(f['downstream']/DAY, 1)
        self.assert_balanced(1e6, 0, end, f)

    def test_transfer_loss_is_explicit_and_water_limited(self):
        d = Design(seepage_fraction_per_day=0, delivery_efficiency=.8)
        end, f = step(100, 0, 0, 0, False, d)
        self.assertEqual(end, 0)
        self.assertEqual(f['release'], 100)
        self.assertEqual(f['delivered'], 80)
        self.assertEqual(f['transfer_loss'], 20)
        self.assert_balanced(100, 0, end, f)

    def test_rainfall_spill_and_evaporation_are_budgeted(self):
        d = Design(capacity_Mm3=1, seepage_fraction_per_day=0)
        end, f = step(1e6, 1, 100, 0, True, d)
        self.assertEqual(end, 1e6)
        self.assertEqual(f['rain'], 10000)
        self.assertEqual(f['spill'], 10000)
        self.assert_balanced(1e6, 1, end, f)
        end, f = step(1, 0, 0, 1e9, False, d)
        self.assertEqual(f['evaporation'], 1)
        self.assertEqual(end+f['release'], 0)
        self.assert_balanced(1, 0, end, f)

    def test_carryover_is_continuous_across_hydrological_years(self):
        data = pd.DataFrame({'Q_m3s':[0,0], 'rain_mm':[0,0], 'pet_Singer_mm':[0,0]},
                            index=pd.to_datetime(['2000-05-31','2000-06-01']))
        d = Design(initial_fraction=.5, seepage_fraction_per_day=0, delivery_efficiency=1)
        run = simulate(data, d)
        self.assertEqual(run.storage_start_m3.iloc[1], run.storage_end_m3.iloc[0])
        self.assertEqual(run.release_m3.iloc[1], 0)

    def test_random_steps_conserve_mass(self):
        rng = np.random.default_rng(766)
        for _ in range(200):
            d = Design(capacity_Mm3=float(rng.uniform(.1,20)), capture_fraction=float(rng.uniform(0,1)),
                       delivery_efficiency=float(rng.uniform(.5,1)), seepage_fraction_per_day=float(rng.uniform(0,.1)))
            initial = rng.uniform(0,d.capacity_Mm3*1e6)
            q = rng.uniform(0,200)
            end, f = step(initial,q,rng.uniform(0,100),rng.uniform(0,15),bool(rng.integers(2)),d)
            self.assert_balanced(initial,q,end,f)
            self.assertLessEqual(end,d.capacity_Mm3*1e6)


if __name__ == '__main__':
    unittest.main()
