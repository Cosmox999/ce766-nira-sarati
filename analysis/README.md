# Historical analysis notes

Author: Ganesh Sonawane

## Scope and inputs

The study covers 1 June 1981-31 May 2011 at Sarati, Nira River. Hydrological years use their starting year: HY 2003 means June 2003-May 2004. All 10,957 daily pairs are present in this window. Missing discharge elsewhere remains missing. The project uses observed discharge, without the CAMELS-IND LSTM predictions.

`preliminary_analysis.py` and `results/` reproduce the preliminary submission. `final_analysis.py` reads the same source extracts and writes `final_results/`. The two verification scripts use separate raw-CSV calculations. The tests in `tests/` check statistical and date-boundary behavior.

## Rainfall standardization

```text
R_k(t) = sum_(i=0)^(k-1) P_(t-i)
SPI_k(t) = Phi_inverse[q + (1-q) G(R_k(t); alpha, beta)]
```

P is monthly rainfall in mm; k is 3 or 6. A separate zero-location gamma distribution is fitted by maximum likelihood to positive totals for each ending calendar month. Alpha is shape, beta is scale in mm, and q is the fraction of zero reference totals. The 1981-2020 reference gives 40 observations for each of 24 fits. Rainfall from 1980 supplies antecedent months. Probabilities are clipped at 1e-10 and 1-1e-10 before the inverse normal transformation.

This follows the probability-standardization approach in the [WMO SPI guide](https://www.droughtmanagement.info/literature/WMO_standardized_precipitation_index_user_guide_en_2012.pdf). Gamma maximum likelihood is the specific distribution/estimator chosen here. Normally dry windows need caution.

The 1981-2020 reference includes years after the discharge-study window. These are retrospective standardized anomalies, as specified in the proposal. An operational forecast evaluation would need a reference fixed before the evaluation period and separate test years.

The source rainfall has at most two decimal places. Monthly and rolling totals retain that 0.01 mm precision before transformation, preserving ties. SPI and discharge are not rounded for calculation. This corrects three preliminary exploratory coefficients: SPI-3 February at zero lead, March at one-month lead, and May at three-month lead. The largest change is about 0.0081. The submitted preliminary statistics, displayed SPI values and September coefficients are unchanged.

## Drought episodes

Episodes are continuous negative SPI runs reaching -1 or below. They end at nonnegative SPI; an exact zero ends the negative run. Duration includes negative lead-in months. Magnitude is minus the sum of SPI during the run. Weaker negative sequences are excluded.

Events are summarized within the study window. Boundary events are flagged as censored if the adjacent SPI remains negative or is unavailable. Their duration and magnitude describe only the observed part. One event per scale continues beyond May 2011. The largest SPI-6 magnitude is a complete event, May 2002-April 2004.

Sensitivity calculations use a 1981-2010 reference and moment-fitted gamma parameters. Event boundaries and counts can change. Within-calendar-month rank correlations are unchanged because these transformations preserve rainfall ranks.

## Recorded flow measures

```text
Z_t(c) = 100 * N_(Q<=c,t) / N_(valid,t)
V_s = 86400 * sum_(d in s) Q_d / 10^6
```

Z is the percentage of valid daily discharge observations at or below c. The main threshold is zero m3/s; sensitivity uses 0.1 m3/s. Q_d is daily mean discharge in m3/s; V_s is volume in million m3 over days in season s. Seasons are June-September, October-November, December-February and March-May.

Preliminary seasonal volumes summed all 30 occurrences. Final seasonal volumes average the 30 occurrences. Dividing the preliminary value by 30 gives the final mean. This is a presentation change.

Positive or missing flow breaks a recorded zero-flow spell. Full-window spells may cross June-May boundaries. Annual spell lengths are truncated at each hydrological-year boundary and labelled accordingly. Monthly fractions pool valid daily observations, rather than averaging percentages with unequal denominators.

## Month-specific associations

```text
rho_(k,l,m) = Corr(rank(SPI_k(t-l)), rank(Z_t(0)))
```

Corr is Pearson correlation of average tied ranks. The rainfall lead l is 0, 1, 2 or 3 months; m identifies the discharge calendar month. Each coefficient uses 30 paired years. Leads are applied to the continuous monthly series before selecting or resampling years. Using flow percentages instead of fractions leaves ranks unchanged.

April's eight coefficients are undefined because the flow fraction is constant. Undefined tests enter the multiplicity family as p=1. They are displayed as NA, not zero. All 96 planned tests remain in the family.

### Intervals

The primary bootstrap samples 30 complete paired hydrological years with replacement, 10,000 times. Each year retains its calendar months and already-aligned lagged rainfall. The same sampled years are used across cells. Ranks are recalculated within every resample, including repeated years and ties. Limits are the 2.5th and 97.5th percentiles.

Dependence checks sample circular runs of two or three adjacent years until 30 years are obtained. Each length uses 10,000 draws. The circular convention joins the end of the record to the beginning. These checks do not establish stationarity or absence of dependence beyond three years.

The tables record the fraction of defined bootstrap coefficients. Intervals are withheld when fewer than 95% of draws are defined. Sparse dry-season variation makes some resamples constant, affecting May and some longer-block calculations. April has no defined coefficient or interval.

Descriptive monthly and seasonal intervals also resample complete years. Seasonal volume intervals describe the mean of 30 occurrences. These are conditional sampling intervals, without original gauge measurement-error bounds.

### Permutation tests

The primary two-sided test permutes rainfall-year labels relative to flow-year labels 19,999 times. The same permutations are used across cells. The p-value is `(1 + count(|rho_null| >= |rho_observed|)) / (B + 1)`. The added one prevents zero Monte Carlo p-values. Holm's step-down correction covers all 96 tests.

The second test permutes ten nonoverlapping three-year rainfall blocks, retaining order within each block. It uses 19,999 draws and a separate Holm correction over the same family. It assumes exchangeable three-year blocks under the null; the primary test assumes exchangeable single years. Borderline adjusted p-values depend on Monte Carlo error and these assumptions. Neither test establishes causation or forecasting skill.

[The SciPy Spearman documentation](https://docs.scipy.org/doc/scipy-1.13.1/reference/generated/scipy.stats.spearmanr.html) recommends permutation tests for small samples. Ties follow the average-rank convention in [rankdata](https://docs.scipy.org/doc/scipy-1.13.1/reference/generated/scipy.stats.rankdata.html).

## Gamma diagnostics and reference uncertainty

For each of 24 fits, the Kolmogorov-Smirnov distance compares positive reference totals with the fitted gamma CDF. I simulate 999 positive samples of the same size, refit each, and recompute the distance. Refitting accounts for parameter estimation in the null samples. P-values use the added-one rule and are supplied before and after Holm correction across 24 checks. Zero probability is recorded separately; this KS check concerns the positive part. Non-rejection does not validate extreme tails.

For SPI intervals, the 40 reference years are resampled 2,000 times, sharing year indices across ending months and scales. Every sample re-estimates zero frequency and gamma parameters, then transforms fixed study rainfall totals. Percentile limits quantify reference-sample uncertainty conditional on the rainfall product and exchangeable reference years. They exclude measurement error, alternative distributions and all forms of serial dependence. Parameter uncertainty does not change within-month rainfall ranks for Spearman correlation.

## Reproducibility and limits

Seed 7662026 supplies separate child random streams for association and gamma calculations. `requirements.txt` pins analysis packages; `run_config.json` records versions and repetitions actually used. `summary.json` stores input SHA-256 checksums. Numerical tables retain full precision; rounding is for display.

`verify_final_results.py` derives all 96 coefficients from raw rainfall totals using Decimal arithmetic, independently of fitted SPI. It also checks source hashes, all study SPI values through the incomplete-gamma function and a separate normal inverse, seasonal volumes, event magnitudes and interval availability.

This historical analysis provides the rainfall/flow context. The later storage-intervention extension is documented in [water_budget_methods.md](water_budget_methods.md), with explicit action-to-outcome equations and daily water balances. The historical correlations alone do not establish intervention benefits. Recorded zeros still require independent channel checks.
