"""Independent raw CSV verification of report statistics and SPI transformation."""
import csv, datetime as dt, json
from decimal import Decimal
from pathlib import Path
from statistics import NormalDist
import math
from scipy.special import gammainc, digamma

ROOT=Path(__file__).resolve().parents[1]
start=dt.date(1981,6,1);end=dt.date(2011,5,31)
raw=ROOT/'data/raw/sarati_source'
results=ROOT/'analysis/results'
rain={};flow={};monthly_rain={}
with (raw/'catchment_mean_forcings/04049.csv').open(newline='') as f:
    for r in csv.DictReader(f):
        day=dt.date(int(r['year']),int(r['month']),int(r['day']))
        value=Decimal(r['prcp(mm/day)'])
        monthly_rain[(day.year,day.month)]=monthly_rain.get((day.year,day.month),Decimal(0))+value
        if start<=day<=end:rain[day]=value
with (raw/'streamflow_timeseries/streamflow_observed.csv').open(newline='') as f:
    for r in csv.DictReader(f):
        day=dt.date(int(r['year']),int(r['month']),int(r['day']))
        if start<=day<=end:flow[day]=Decimal(r['4049'])
summary=json.loads((results/'summary.json').read_text())
assert len(rain)==len(flow)==(end-start).days+1==summary['n_days']
assert all(v.is_finite() and v>=0 for v in rain.values())
assert all(v.is_finite() and v>=0 for v in flow.values())
assert sum(v==0 for v in flow.values())==summary['zero_days']==6984
assert sum(v<=Decimal('.1') for v in flow.values())==summary['nearzero_days']==7016
assert abs(float(sum(rain.values())/30)-summary['mean_annual_rain_mm'])<1e-9
cur=longest=0;run_end=None
for day in sorted(flow):
    cur=cur+1 if flow[day]==0 else 0
    if cur>longest:longest=cur;run_end=day
assert longest==summary['longest_zero_run']['days']==293
assert str(run_end)==summary['longest_zero_run']['end']
assert str(run_end-dt.timedelta(days=longest-1))==summary['longest_zero_run']['start']
with (results/'seasonal_metrics.csv').open(newline='') as f:
    for row in csv.DictReader(f):
        months={'Jun-Sep':[6,7,8,9],'Oct-Nov':[10,11],'Dec-Feb':[12,1,2],'Mar-May':[3,4,5]}[row['season']]
        vals=[v for day,v in flow.items() if day.month in months]
        assert len(vals)==int(row['valid_Q'])
        assert sum(v==0 for v in vals)==int(row['zero_days'])
        assert abs(sum(v==0 for v in vals)/len(vals)-float(row['zero_fraction']))<1e-12
        volume=sum(vals)*Decimal(86400)/Decimal(1000000)
        assert abs(float(volume)-float(row['observed_volume_Mm3']))<1e-6

# Rebuild reported September 2003 SPI windows without pandas rolling or
# scipy.stats gamma/normal distribution wrappers. Check MLE score equation.
fits=list(csv.DictReader((results/'spi_gamma_fits.csv').open()))
monthly=list(csv.DictReader((results/'monthly_spi_and_flow.csv').open()))
target=next(r for r in monthly if r['month']=='2003-09-01')
for k in [3,6]:
    ref=[]
    for year in range(1981,2021):
        ref.append(float(sum(monthly_rain[(year,m)] for m in range(10-k,10))))
    positive=[x for x in ref if x>0]
    fit=next(r for r in fits if int(r['scale_months'])==k and r['ending_month']=='9')
    shape=float(fit['shape']);scale=float(fit['scale'])
    assert abs(scale-sum(positive)/len(positive)/shape)<1e-10
    score=math.log(shape)-digamma(shape)-math.log(sum(positive)/len(positive))+sum(math.log(x) for x in positive)/len(positive)
    assert abs(score)<1e-8
    total=float(sum(monthly_rain[(2003,m)] for m in range(10-k,10)))
    p0=sum(x==0 for x in ref)/len(ref)
    h=p0+(1-p0)*gammainc(shape,total/scale)
    spi=NormalDist().inv_cdf(h)
    assert abs(spi-float(target[f'SPI{k}']))<1e-10
print('PASS: raw daily counts, rainfall totals, seasonal fractions/volumes, longest spell, and reported SPI-3/SPI-6.')
(results/'verification.txt').write_text('PASS: independent raw-CSV checks of counts, rainfall, seasonal metrics and longest spell.\nPASS: independent September 2003 SPI transformation and gamma MLE score checks.\n')
