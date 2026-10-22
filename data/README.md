# Data provenance and attribution

These are station-specific extracts of **CAMELS-IND version 2.2**, published on 13 March 2025 and accessed on 21 October 2026.

- Repository: https://zenodo.org/records/14999580
- Version DOI: https://doi.org/10.5281/zenodo.14999580
- Original archive: `CAMELS_IND_Catchments_Streamflow_Sufficient.zip`
- Original archive size: 178,141,731 bytes
- Archive MD5 published by Zenodo and verified after download: `3993c25ba7d7b86df0541de91e094f39`
- Dataset licence in the Zenodo record: [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).

The dataset is credited to **Nikunj K. Mangukiya, Kanneganti Bhargav Kumar, Pankaj Dey, Shailza Sharma, Vijaykumar Bejagam, Pradeep P. Mujumdar and Ashutosh Sharma**. Original rainfall observations are from IMD; observed discharge was compiled from CWC/India-WRIS by the dataset authors. Ganesh Sonawane is the author of this project, rather than the original hydrometeorological dataset.

## Included source extracts

All paths below are relative to `raw/sarati_source/` and retain the original archive filenames.

| File | Selection | Units |
| --- | --- | --- |
| `catchment_mean_forcings/04049.csv` | Original `year`, `month`, `day`, `prcp(mm/day)` columns | Catchment-mean daily rainfall, mm/day |
| `streamflow_timeseries/streamflow_observed.csv` | Original `year`, `month`, `day`, `4049` columns | Observed daily mean discharge, m3/s |
| `attributes_csv/camels_ind_name.csv` | Original row for gauge 4049 | Gauge identity and river information |
| `attributes_csv/camels_ind_topo.csv` | Original row for gauge 4049 | Includes coordinates and catchment area |

Daily files cover 1 January 1980-31 December 2020, including missing discharge outside the study window. Original numeric strings and empty cells were retained; unrelated columns and stations were omitted. CSV formatting was normalized to UTF-8 with LF line endings. No rainfall or discharge values were imputed, rounded again or replaced by LSTM estimates.

The delineated GHI catchment area is 6,807.1 km2; the metadata also contain a separate reported CWC area of 7,200 km2. This project uses the delineated area consistent with its proposal. Station 04049 in filenames is represented as 4049 in table columns.

`analysis/results/input_sha256.json` records hashes of the extracts used in the reproduced calculation. These hashes describe the selected files in this repository, rather than the larger original multi-station files.

## Citation

Mangukiya, N. K., et al. (2025). *CAMELS-IND: hydrometeorological time series and catchment attributes for 472 catchments in Peninsular India* (version 2.2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14999580

Associated paper: Mangukiya, N. K., et al. (2025). CAMELS-IND: hydrometeorological time series and catchment attributes for 228 catchments in Peninsular India. *Earth System Science Data*, 17, 461-491. https://doi.org/10.5194/essd-17-461-2025
