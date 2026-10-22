# Data sources

The files in this folder are Sarati extracts from CAMELS-IND version 2.2. This version was published on 13 March 2025. I downloaded the source data on 21 October 2026.

- Repository: https://zenodo.org/records/14999580
- Version DOI: https://doi.org/10.5281/zenodo.14999580
- Original archive: `CAMELS_IND_Catchments_Streamflow_Sufficient.zip`
- Original archive size: 178,141,731 bytes
- Archive MD5 published by Zenodo and verified after download: `3993c25ba7d7b86df0541de91e094f39`
- Dataset licence in the Zenodo record: [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).

CAMELS-IND was compiled by Nikunj K. Mangukiya, Kanneganti Bhargav Kumar, Pankaj Dey, Shailza Sharma, Vijaykumar Bejagam, Pradeep P. Mujumdar and Ashutosh Sharma. The rainfall data come from IMD, and the observed discharge data were compiled from CWC/India-WRIS.

## Files used

The paths below are inside `raw/sarati_source/`. I kept the original filenames so they can be matched to the source archive.

| File | Selection | Units |
| --- | --- | --- |
| `catchment_mean_forcings/04049.csv` | Original `year`, `month`, `day`, `prcp(mm/day)` columns | Catchment-mean daily rainfall, mm/day |
| `streamflow_timeseries/streamflow_observed.csv` | Original `year`, `month`, `day`, `4049` columns | Observed daily mean discharge, m3/s |
| `attributes_csv/camels_ind_name.csv` | Original row for gauge 4049 | Gauge identity and river information |
| `attributes_csv/camels_ind_topo.csv` | Original row for gauge 4049 | Includes coordinates and catchment area |

The daily files cover 1 January 1980-31 December 2020. I selected the required columns and Sarati metadata rows and saved them as UTF-8 CSV files with LF line endings. The original numeric values and empty cells were kept. Missing discharge outside the study period was left missing, and the analysis uses observed discharge throughout.

There are two catchment areas in the metadata: the delineated GHI area is 6,807.1 km2, while the reported CWC area is 7,200 km2. I use the delineated area, as in the proposal. Station 04049 in the filenames appears as 4049 in the table columns.

`analysis/results/input_sha256.json` contains the checksums for the extracts in this repository. The MD5 listed above is for the complete archive downloaded from Zenodo.

## Citation

Mangukiya, N. K., et al. (2025). *CAMELS-IND: hydrometeorological time series and catchment attributes for 472 catchments in Peninsular India* (version 2.2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14999580

Associated paper: Mangukiya, N. K., et al. (2025). CAMELS-IND: hydrometeorological time series and catchment attributes for 228 catchments in Peninsular India. *Earth System Science Data*, 17, 461-491. https://doi.org/10.5194/essd-17-461-2025
