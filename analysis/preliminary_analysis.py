"""Calculate rainfall, SPI and recorded zero-flow summaries for Sarati."""
# Author: Ganesh Sonawane
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw/sarati_source'
OUT = ROOT / 'analysis/results'
OUT.mkdir(parents=True, exist_ok=True)
P = pd.read_csv(RAW / 'catchment_mean_forcings/04049.csv')
Q = pd.read_csv(RAW / 'streamflow_timeseries/streamflow_observed.csv', usecols=['year','month','day','4049'])
for df in [P,Q]:
    # Check the full daily date sequence before selecting the study period.
    df.index = pd.to_datetime(df[['year','month','day']])
    assert not df.index.duplicated().any()
    assert df.index.equals(pd.date_range('1980-01-01','2020-12-31'))
d = pd.DataFrame({'P_mm':P['prcp(mm/day)'], 'Q_m3s':Q['4049']})
assert not (d[['P_mm','Q_m3s']] < 0).any().any(), 'Investigate negative/sentinel values'
assert not np.isinf(d.to_numpy()).any()
a = d.loc['1981-06-01':'2011-05-31'].copy()
# Label each June-May hydrological year by its starting year.
a['HY'] = a.index.year - (a.index.month < 6)
a['season'] = np.select([a.index.month.isin([6,7,8,9]),a.index.month.isin([10,11]),a.index.month.isin([12,1,2])],['Jun-Sep','Oct-Nov','Dec-Feb'],default='Mar-May')
assert len(a)==10957
def longest(x):
    # A positive value or missing observation ends a recorded zero-flow spell.
    flag = x.eq(0)
    return int(flag.groupby((~flag).cumsum()).sum().max())

annual=[]
for y,g in a.groupby('HY'):
    annual.append(dict(HY=int(y), days=len(g), valid_Q=int(g.Q_m3s.count()), missing_Q=int(g.Q_m3s.isna().sum()), rainfall_mm=g.P_mm.sum(min_count=len(g)), zero_days=int(g.Q_m3s.eq(0).sum()), nearzero_days=int(g.Q_m3s.le(.1).sum()), longest_recorded_zero_days=longest(g.Q_m3s)))
annual=pd.DataFrame(annual)
seasonal=[]
for s in ['Jun-Sep','Oct-Nov','Dec-Feb','Mar-May']:
    g=a[a.season.eq(s)]
    n=int(g.Q_m3s.count())
    seasonal.append(dict(season=s,calendar_days=len(g),valid_Q=n,missing_Q=len(g)-n,mean_rainfall_mm=g.P_mm.sum()/30,zero_days=int(g.Q_m3s.eq(0).sum()),zero_fraction=g.Q_m3s.eq(0).sum()/n,nearzero_fraction=g.Q_m3s.le(.1).sum()/n,observed_volume_Mm3=g.Q_m3s.sum()*86400/1e6))
seasonal=pd.DataFrame(seasonal)
monthly=d.resample('MS').agg(P_mm=('P_mm','sum'),P_n=('P_mm','count'),Q_n=('Q_m3s','count'),zero_n=('Q_m3s',lambda x:x.eq(0).sum()))
monthly.loc[monthly.P_n.ne(monthly.index.days_in_month),'P_mm']=np.nan
monthly['zero_fraction']=monthly.zero_n/monthly.Q_n.replace(0,np.nan)
fit_records=[]
for k in [3,6]:
    roll=monthly.P_mm.rolling(k,min_periods=k).sum()
    spi=pd.Series(np.nan,index=monthly.index)
    for m in range(1,13):
        # Compare rainfall totals ending in the same calendar month.
        ref=roll.loc['1981-01-01':'2020-12-01']
        ref=ref[ref.index.month==m].dropna()
        pos=ref[ref>0]
        assert len(ref)==40 and len(pos)>10
        shape,loc,scale=stats.gamma.fit(pos,floc=0)
        # Include zero totals separately from the gamma fit to positive totals.
        p0=float(ref.eq(0).mean())
        use=roll[roll.index.month==m].dropna()
        h=p0+(1-p0)*stats.gamma.cdf(use,shape,loc=0,scale=scale)
        spi.loc[use.index]=stats.norm.ppf(np.clip(h,1e-10,1-1e-10))
        fit_records.append(dict(scale_months=k,ending_month=m,n_reference=len(ref),zero_probability=p0,shape=shape,scale=scale))
    monthly[f'SPI{k}']=spi

# Calculate correlations separately for each month to account for seasonality.
# A positive lag pairs earlier rainfall SPI with the current flow fraction.
correlations=[]
for k in [3,6]:
    for lag in range(4):
        pair=monthly[[f'SPI{k}','zero_fraction','Q_n']].copy()
        pair['spi_lead']=pair[f'SPI{k}'].shift(lag)
        pair=pair.loc['1981-06-01':'2011-05-01']
        # Keep only months with discharge recorded on every day.
        pair=pair[pair.Q_n.eq(pair.index.days_in_month)]
        for m in range(1,13):
            g=pair[pair.index.month==m][['spi_lead','zero_fraction']].dropna()
            rho=np.nan
            if len(g)>=10 and g.spi_lead.nunique()>1 and g.zero_fraction.nunique()>1:
                rho=float(stats.spearmanr(g.spi_lead,g.zero_fraction).statistic)
            correlations.append(dict(spi_scale=k,rainfall_lead_months=lag,flow_month=m,n=len(g),spearman_rho=rho))

flag=a.Q_m3s.eq(0)
# Allow these spells to continue across hydrological-year boundaries.
runs=[]
for _,g in a.groupby((~flag).cumsum()):
    z=g[g.Q_m3s.eq(0)]
    if len(z):runs.append(dict(start=str(z.index[0].date()),end=str(z.index[-1].date()),days=len(z)))
runs=sorted(runs,key=lambda x:x['days'],reverse=True)
summary=dict(n_days=len(a),rainfall_missing=int(a.P_mm.isna().sum()),flow_missing=int(a.Q_m3s.isna().sum()),flow_valid=int(a.Q_m3s.count()),zero_days=int(a.Q_m3s.eq(0).sum()),nearzero_days=int(a.Q_m3s.le(.1).sum()),zero_fraction=float(a.Q_m3s.eq(0).sum()/a.Q_m3s.count()),nearzero_fraction=float(a.Q_m3s.le(.1).sum()/a.Q_m3s.count()),mean_annual_rain_mm=float(annual.rainfall_mm.mean()),min_rain_HY=int(annual.loc[annual.rainfall_mm.idxmin(),'HY']),min_rain_mm=float(annual.rainfall_mm.min()),max_rain_HY=int(annual.loc[annual.rainfall_mm.idxmax(),'HY']),max_rain_mm=float(annual.rainfall_mm.max()),longest_zero_run=runs[0],complete_flow_HYs=int(annual.missing_Q.eq(0).sum()),full_record_flow_valid=int(d.Q_m3s.count()),full_record_days=len(d),analysis_first=str(a.index.min().date()),analysis_last=str(a.index.max().date()))
annual.to_csv(OUT/'annual_metrics.csv',index=False)
seasonal.to_csv(OUT/'seasonal_metrics.csv',index=False)
monthly.to_csv(OUT/'monthly_spi_and_flow.csv',index_label='month')
pd.DataFrame(fit_records).to_csv(OUT/'spi_gamma_fits.csv',index=False)
pd.DataFrame(correlations).to_csv(OUT/'exploratory_correlations.csv',index=False)
pd.DataFrame(runs).to_csv(OUT/'zero_flow_runs.csv',index=False)
a.to_csv(OUT/'sarati_analysis_daily.csv',index_label='date')
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
hashes={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in RAW.rglob('*') if p.is_file()}
(OUT/'input_sha256.json').write_text(json.dumps(hashes,indent=2))

# Show the monthly pattern with June at the start of the axis.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(7.0,1.8),layout='constrained')
months=[6,7,8,9,10,11,12,1,2,3,4,5]
rain=a.groupby(a.index.month).P_mm.sum()/30
frac=a.groupby(a.index.month).apply(lambda g:g.Q_m3s.eq(0).sum()/g.Q_m3s.count())*100
labels=['Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May']
axs[0].bar(range(12),rain.reindex(months),color='#2C7A9B',width=.75)
axs[1].bar(range(12),frac.reindex(months),color='#BA6A35',width=.75)
for ax in axs:
    ax.set_xticks(range(12),labels)
    ax.grid(axis='y',alpha=.18)
    ax.set_axisbelow(True)
    ax.tick_params(axis='x',labelsize=7)
    ax.tick_params(axis='y',labelsize=8)
axs[0].set_title('(a) Mean monthly rainfall',loc='left',fontsize=10)
axs[0].set_ylabel('mm/month',fontsize=9)
axs[1].set_title('(b) Recorded zero-flow days',loc='left',fontsize=10)
axs[1].set_ylabel('% of valid observations',fontsize=9)
axs[1].set_ylim(0,100)
fig.savefig(OUT/'preliminary_figure.png',dpi=250)
fig.savefig(OUT/'preliminary_figure.pdf')
plt.close(fig)
print(json.dumps(summary,indent=2))
print(seasonal.to_string(index=False))
print(annual.to_string(index=False))
