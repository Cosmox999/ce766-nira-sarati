# Rainfall deficits and river-flow intermittency at Sarati

Author: Ganesh Sonawane

This is my term project for CE 766, Watershed Management. I am studying the Nira River catchment upstream of Sarati in Maharashtra to see how rainfall deficits are associated with periods of recorded zero flow.

The study covers June 1981-May 2011, giving 30 complete hydrological years. I am using CAMELS-IND station 04049, which appears as 4049 in the discharge file. The delineated catchment area is 6,807.1 km2.

## Work completed so far

I have checked the daily rainfall and discharge records, calculated monthly and seasonal summaries, identified continuous zero-flow spells, and calculated SPI-3 and SPI-6. I have also made an initial set of rainfall-flow correlation calculations, but the uncertainty analysis is still pending.

Both daily series are complete within the study period. The main results so far are:

| Measure | Result |
| --- | --- |
| Days covered | 10,957 |
| Days with recorded zero discharge | 6,984 (63.74%) |
| Days with discharge at or below 0.1 m3/s | 7,016 (64.03%) |
| Zero-flow fraction during March-May | 98.84% |
| Mean rainfall per hydrological year | 736.85 mm |
| Share of rainfall received during June-September | 79.26% |
| Longest recorded zero-flow spell | 293 days, from 23 October 2002 to 11 August 2003 |

![Monthly rainfall and recorded zero-flow frequency at Sarati](analysis/results/preliminary_figure.png)

The rainfall bars are monthly totals averaged over the study period. The zero-flow bars show the percentage of valid daily discharge observations recorded as zero in each calendar month. The months run from June to May.

## Data

The [data folder](data/raw/sarati_source) contains the rainfall, observed discharge and station details used in the analysis. These are Sarati extracts from [CAMELS-IND v2.2](https://zenodo.org/records/14999580), taken from `CAMELS_IND_Catchments_Streamflow_Sufficient.zip`.

I kept the original values and missing entries and selected the columns needed for this project. Missing discharge outside the study period remains missing. The [data notes](data/README.md) list the source files, units and attribution.

## Running the scripts

Use Python 3.11 or 3.12. From the repository folder, create a virtual environment:

```bash
python -m venv .venv
```

Activate it in PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

On Linux or macOS, use:

```bash
source .venv/bin/activate
```

Install the packages and run the analysis, followed by the checks:

```bash
python -m pip install -r requirements.txt
python analysis/preliminary_analysis.py
python analysis/verify_results.py
```

The first script writes the tables and figure to [analysis/results](analysis/results). The second checks the reported totals, seasonal measures, longest zero-flow spell and the September 2003 SPI values using a separate calculation. These commands also run through GitHub Actions when changes are pushed.

## A few points about the method

Hydrological years run from June to May. For example, HY 2003 means June 2003-May 2004. Zero-flow fractions use the number of valid discharge observations as the denominator.

SPI-3 and SPI-6 use three- and six-month rainfall totals. The gamma distributions are fitted separately for each ending calendar month, using 1981-2020 as the reference period. Rainfall from 1980 supplies the preceding months needed at the start of the reference period.

Rainfall-flow correlations are calculated separately for each calendar month, with rainfall leads of zero to three months. This avoids combining the usual seasonal cycle with year-to-year differences. More details are in the [analysis notes](analysis/README.md).

## Next steps

I still need to examine the reliability of the recorded zeros, check the SPI fits, summarize rainfall-deficit episodes and estimate uncertainty using hydrological-year block bootstrapping. Reservoir operation and abstraction may also affect discharge at Sarati. These results therefore describe the available records; they do not establish the cause of the zero-flow conditions or the benefit of a particular intervention.

## References

1. Mangukiya, N. K., et al. (2025). *CAMELS-IND: hydrometeorological time series and catchment attributes for 472 catchments in Peninsular India*, version 2.2. Zenodo. https://doi.org/10.5281/zenodo.14999580
2. Mangukiya, N. K., et al. (2025). CAMELS-IND: hydrometeorological time series and catchment attributes for 228 catchments in Peninsular India. *Earth System Science Data*, 17, 461-491. https://doi.org/10.5194/essd-17-461-2025
3. World Meteorological Organization (2012). *Standardized Precipitation Index User Guide* (M. Svoboda, M. Hayes and D. Wood), WMO-No. 1090. [Guide](https://www.droughtmanagement.info/literature/WMO_standardized_precipitation_index_user_guide_en_2012.pdf).
