# Analysis notes

`preliminary_analysis.py` calculates the preliminary results from the Sarati data files. `verify_results.py` checks the main numbers using a separate calculation.

## Reading the tables

The study period is 1 June 1981-31 May 2011, including both dates. There are 10,957 days, with no missing rainfall or discharge observations in this period. The longer discharge file does have missing entries outside this window.

`HY` refers to the year in which the hydrological year starts. HY 2003, for example, covers June 2003-May 2004. Leap days are included.

The zero-flow fraction is the number of valid daily discharge values equal to zero, divided by the number of valid observations. I also repeated the count using a threshold of 0.1 m3/s to check how much small positive values affect the result.

The seasonal groups are June-September, October-November, December-February and March-May. Rainfall in the seasonal table is averaged over 30 occurrences of each season. Discharge volume is summed over all 30 occurrences and reported in million m3. Divide that volume by 30 to obtain the average volume for one season. The conversion is:

```text
volume in million m3 = sum of daily mean discharge in m3/s x 86,400 / 1,000,000
```

The annual table ends each zero-flow spell at the June-May year boundary. The full spell table allows a spell to continue across that boundary. This is why the 293-day longest spell in the full record is longer than the maxima within individual hydrological years.

## Calculating SPI

SPI-k uses the rainfall total for k months, including the ending month. The reference period consists of ending months in 1981-2020, giving 40 samples for each calendar month. The preceding months needed for early 1981 come from the 1980 rainfall data.

For each ending month, a gamma distribution is fitted to positive rainfall totals by maximum likelihood, with its location fixed at zero. The probability of a zero total is included separately. If p0 is that probability and G(x) is the fitted gamma CDF, the combined probability is:

```text
H(x) = p0 + (1 - p0) G(x)
```

SPI is obtained by applying the inverse standard-normal CDF to H(x). Probabilities are limited to the range 1e-10 to 1-1e-10 to avoid infinite values. This limit does not affect the September 2003 values reported in the preliminary submission.

The fits still need further checking, especially for normally dry months. The verification script checks the gamma fit equations and the conversion to SPI for September 2003 at both accumulation periods.

## Rainfall-flow comparisons

The correlation file contains 96 initial comparisons: two SPI periods, four rainfall leads and twelve calendar months. A lead of l months pairs SPI(t-l) with the zero-flow fraction in month t. Only months with complete discharge observations are used.

There are at most 30 pairs for each calendar-month correlation. Where a series is constant, the correlation is left blank because it is undefined. Confidence intervals, dependence between years and the effect of testing several comparisons still need to be assessed before these coefficients are interpreted.

## Checking the results

The verification script reads the CSV files with Python's CSV reader and uses Decimal arithmetic for the daily totals. It checks the observation counts, rainfall totals, seasonal fractions and volumes, longest spell and its dates, and the reported SPI calculations. A successful run writes `results/verification.txt`.

These checks confirm the calculations against the supplied files. Checking the original gauge records and reservoir operations is a separate part of the remaining work.
