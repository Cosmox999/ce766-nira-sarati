"""Verify exported balances and reconstruct the central run independently."""
# Author: Ganesh Sonawane
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'analysis/water_budget_results'


def main():
    summary = json.loads((OUT/'summary.json').read_text())
    for name, value in summary['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == value
    for path in OUT.glob('daily_*.csv'):
        with path.open(newline='') as file:
            previous = None
            for row in csv.DictReader(file):
                value = {key:Decimal(text) for key,text in row.items() if key!='date'}
                delta = value['storage_end_m3']-value['storage_start_m3']
                inputs = value['diversion_m3']+value['rain_m3']
                outputs = sum(value[key] for key in ['evaporation_m3','seepage_m3','release_m3','spill_m3'])
                assert abs(inputs-outputs-delta) < Decimal('0.00001'),(path,row['date'])
                net = value['inflow_m3']+value['rain_m3']-sum(value[key] for key in
                      ['downstream_m3','evaporation_m3','seepage_m3','transfer_loss_m3'])
                assert abs(net-delta) < Decimal('0.00001')
                if previous is not None:
                    assert previous == value['storage_start_m3']
                previous = value['storage_end_m3']
    weather = pd.read_csv(ROOT/'data/raw/water_budget_forcings/04049.csv')
    previous_rain = pd.read_csv(ROOT/'data/raw/sarati_source/catchment_mean_forcings/04049.csv')
    np.testing.assert_array_equal(weather['prcp(mm/day)'], previous_rain['prcp(mm/day)'])
    flow = pd.read_csv(ROOT/'data/raw/sarati_source/streamflow_timeseries/streamflow_observed.csv')
    weather.index = flow.index = pd.to_datetime(weather[['year','month','day']])
    model = pd.read_csv(OUT/'daily_K10_f10.csv',index_col=0,parse_dates=True)
    storage = 0.0
    for date,row in model.iterrows():
        q = float(flow.loc[date,'4049'])
        summer = date.month in [6,7,8,9]
        capture = min(max(q-1,0)*.1*86400,5*86400,10e6-storage) if summer else 0.0
        volume = storage+capture
        area = 1e6*(volume/10e6)**(2/3)
        rain = float(weather.loc[date,'prcp(mm/day)'])*area/1000
        evap = min(float(weather.loc[date,'pet(mm/day)'])*area/1000,volume+rain)
        remaining = volume+rain-evap
        seepage = remaining*.001
        remaining -= seepage
        release = min(max(1-q,0)*86400/.9,2*86400,remaining) if not summer else 0.0
        remaining -= release
        spill = max(remaining-10e6,0)
        storage = remaining-spill
        downstream = q*86400-capture+.9*release+spill
        for name,expected in [('storage_end_m3',storage),('diversion_m3',capture),
                              ('evaporation_m3',evap),('seepage_m3',seepage),
                              ('release_m3',release),('downstream_m3',downstream)]:
            assert np.isclose(row[name],expected,rtol=1e-11,atol=1e-5),(date,name,row[name],expected)
    nominal = pd.read_csv(OUT/'nominal_designs.csv')
    ensemble = pd.read_csv(OUT/'uncertainty_ensemble.csv')
    assert len(nominal)==9 and len(ensemble)==324
    assert ensemble.max_system_residual_m3.max()<1e-5
    assert ensemble.max_reservoir_residual_m3.max()<1e-5
    baseline = pd.read_csv(OUT/'daily_baseline.csv')
    np.testing.assert_allclose(baseline.inflow_m3,baseline.downstream_m3,rtol=0,atol=0)
    rows = pd.read_csv(OUT/'annual_storage_ledgers.csv')
    assert len(rows)==270
    for _,r in rows.iterrows():
        balance = r.diversion_Mm3+r.rain_Mm3-r.evaporation_Mm3-r.seepage_Mm3-r.release_Mm3-r.spill_Mm3
        assert abs(balance-(r.storage_end_Mm3-r.storage_start_Mm3))<1e-8
    result = ('PASS: source hashes; daily reservoir and combined-system closure for baseline and all nine nominal runs; '
              'continuous carryover; independent reconstruction of every central-run day; '
              '324 ensemble constraints; 270 annual budgets.\n')
    (OUT/'verification.txt').write_text(result,encoding='utf-8')
    print(result.strip())


if __name__ == '__main__':
    main()
