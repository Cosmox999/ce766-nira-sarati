# Analysis definitions and verification

`preliminary_analysis.py` reads the included Sarati data extracts and recreates the preliminary results. `verify_results.py` checks the principal statistics through a separate calculation.

## Definitions

- Study window: 1 June 1981-31 May 2011 inclusive, comprising 10,957 days. Both input series are complete in this window. This completeness statement does not apply to the entire 1980-2020 discharge file.
- `HY` is the starting year of a June-May hydrological year. HY 2003 means June 2003-May 2004. Leap days are retained.
- Recorded zero-flow fraction = valid daily discharge observations equal to zero / valid daily discharge observations. The near-zero sensitivity threshold is 0.1 m3/s. Missing observations are not counted as zero.
- Seasons: June-September, October-November, December-February and March-May. Seasonal rainfall is averaged over 30 occurrences. Seasonal volume is the sum over all 30 occurrences, in million m3; divide by 30 for the average volume per season.
- Volume in million m3 = sum of daily mean discharge in m3/s x 86,400 / 1,000,000.
- The figure averages monthly rainfall totals across the study window and pools valid daily zero-flow observations by calendar month. It does not average monthly percentages with unequal denominators.
- Annual zero-flow spell lengths are truncated at June-May year boundaries. The full-window spell table permits a spell to cross a hydrological-year boundary. Consequently the reported 293-day full-window maximum differs from the annual truncated maxima.

## Standardized Precipitation Index

SPI-k uses the rainfall sum for k months including the ending month. Reference samples have ending months in 1981-2020, giving 40 samples for each ending calendar month; 1980 rainfall supplies the required antecedents.

Fit gamma shape and scale by maximum likelihood to positive totals, with gamma location fixed at zero. Let p0 be the reference probability of a zero total and G(x) the fitted gamma CDF. Then H(x) = p0 + (1-p0)G(x), and SPI is the inverse standard-normal CDF of H(x). CDF clipping at 1e-10 and 1-1e-10 prevents infinite tails. It does not affect the reported September 2003 values.

Gamma-fit adequacy and sensitivity in normally dry months remain provisional. The independent verification checks the September 2003 SPI-3 and SPI-6 transformation and the gamma maximum-likelihood score equations.

## Exploratory correlations

`exploratory_correlations.csv` contains 96 screening comparisons: two SPI scales x four rainfall leads x twelve calendar months. A positive lead l pairs SPI(t-l) with zero-flow fraction(t). Only complete observed discharge months enter the calculation. A constant series has an undefined correlation, recorded as blank.

These coefficients have not yet received block-bootstrap intervals, dependence assessment or a multiple-comparison assessment. There are at most 30 paired observations per calendar-month coefficient in the study period. They must not be interpreted as established relationships or causal effects.

## Verification

The verification script uses the standard-library CSV reader and Decimal arithmetic independently of the pandas aggregations. It checks daily counts, rainfall totals, seasonal fractions and volumes, the longest spell and its endpoints, and the reported SPI transformations. A successful run writes `results/verification.txt`.

The checks establish reproducibility of these calculations. They do not independently validate the original gauge recording practices, upstream operations or the underlying rainfall product.
