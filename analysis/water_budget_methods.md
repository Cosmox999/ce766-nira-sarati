# Seasonal storage water budget

Author: Ganesh Sonawane

## System boundary

This is a hypothetical intervention-screening model for the Nira catchment at Sarati. The historical rainfall/SPI and flow analysis remains the baseline context. A new off-channel storage unit is placed conceptually just downstream of the gauge. Observed regulated flow is the fixed upstream boundary. It is not natural inflow to Bhatghar or Vir, and the model does not reconstruct or reoperate those reservoirs.

A diversion junction splits the measured daily inflow into bypass flow and storage inflow. Storage returns water to the river during October-May. Seepage and delivery losses leave the modeled reach-storage system; no same-day groundwater return is credited. Rain on the wetted pool is an additional external input. Catchment rainfall runoff is not added to measured flow a second time.

## Daily balances

Volumes are in m3 and flow is in m3/s. S is storage at the start of a step. K is capacity; D diversion; H direct pool rainfall; E evaporation; G seepage; R release from storage; W spill; eta delivery efficiency. T=(1-eta)R is transfer loss.

```text
S_next - S = D + H - E - G - R - W
Q_down = Q_obs - D/86400 + (eta*R + W)/86400
86400*Q_obs + H = 86400*Q_down + E + G + T + (S_next-S)
```

The combined-system equation removes the internal diversion and release transfers. Both equations are checked for every day. Reported closure is numerical bookkeeping, not validation of the assumed hydrology.

The storage starts empty on 1 June 1981 and runs continuously to 31 May 2011. It is never reset at year boundaries. A full-start sensitivity is reported separately.

## Diversion, geometry and release rules

During June-September:

```text
D = min(f * max(Q_obs - q_b, 0) * 86400, q_intake * 86400, K-S)
```

Diversion is zero during the other months. The capture fraction f is 0.05, 0.10 or 0.20. A bypass floor q_b=1 m3/s and intake capacity q_intake=5 m3/s are assumed screening constraints. Existing flows below the floor remain unchanged; no water is diverted from them. The model does not establish that flows above this floor are unallocated or ecologically surplus.

Storage capacity is 5, 10 or 20 million m3. The assumed full-pool mean depth is 10 m, giving full-pool areas of 0.5, 1 and 2 km2. A synthetic area-storage curve is used:

```text
A = (K / 10 m) * ((S+D)/K)^(2/3)
H = P * A / 1000
E = min(c_e * PET * A / 1000, S+D+H)
G = lambda * (S+D+H-E)
```

P and PET are daily depths in mm. Area is evaluated after diversion, once per step. Direct rainfall is limited to the wetted area; runoff from dry reservoir footprint is excluded. Evaporation and seepage cannot exceed available water. These geometry and loss choices are uncalibrated. Seepage is a loss from the surface system, without claiming a usable aquifer benefit.

During October-May:

```text
R = min(max(q_target-Q_obs,0) * 86400 / eta, q_release*86400, S+D+H-E-G)
```

Release is zero during June-September. The assumed downstream target is 1 m3/s, with 0.5 and 2 m3/s sensitivity tests. The outlet limit is 2 m3/s. The target is an illustrative service threshold, not an established ecological requirement. Spill removes any remaining excess over K and returns it to the river. The calculation order is diversion, area, rain, evaporation, seepage, release, spill. Half-day steps hold observed daily means constant and test numerical sensitivity to this ordering.

## Design and uncertainty runs

There are nine capacity-capture combinations. Central assumptions are Singer PET with multiplier 1, lambda=0.001/day and eta=0.9. Every design is evaluated under 36 deterministic combinations: Singer PET times 0.8, 1.0 or 1.2, or GLEAM PET times 1.0; seepage 0.0005, 0.001 or 0.002/day; and delivery efficiency 0.8, 0.9 or 1.0. This gives 324 design-uncertainty combinations. Their ranges are not confidence intervals or probability statements.

For each design, four additional runs use full initial storage, a target of 0.5 m3/s, a target of 2 m3/s, or half-day timesteps. These 36 runs are sensitivity checks. The baseline is pass-through with zero new capacity. The nominal central settings are already included in the 324 combinations.

## Performance measures

October-May deficit volume is the sum of max(q_target-Q_down,0)*86400, averaged over 30 hydrological years and expressed in million m3/year. Deficit reduction is 100*(1-intervention deficit/baseline deficit). Reliability is the fraction of October-May days meeting the target, with a numerical tolerance of 1e-9 m3/s. Diversion is reported both as a volume and as a percentage of observed June-September volume. Capacity is a construction-scale proxy; no economic optimum is claimed.

The 2002-2004 example was identified by the earlier SPI-6 analysis, which found a complete negative-SPI episode during May 2002-April 2004. Its trajectory is shown alongside the same measured rainfall and flow, without assuming a cause for all recorded zeros.

## Catchment-budget context

`catchment_budget_diagnostic.csv` provides P - AET_GLEAM - Q_obs converted to catchment depth over 6,807.1 km2 for each hydrological year. The residual includes unknown storage change, groundwater exchange, transfers, withdrawals, product mismatch and errors. It is not assigned entirely to storage or calibrated away. This separate diagnostic is not the inflow-generation model for the intervention.

## Sources and limitations

Data: [CAMELS-IND v2.2](https://doi.org/10.5281/zenodo.14999580) and its [dataset paper](https://doi.org/10.5194/essd-17-461-2025). Rainfall is based on IMD's 0.25-degree grid. Singer PET is based on a 0.1-degree product; alternative PET/AET are from GLEAM. The repository uses supplied catchment means, not newly extracted grid cells.

The continuity equation and evaporation-volume principle are described in the USACE HEC-HMS Technical Reference Manual: [reservoir concepts](https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/reservoir-modeling/reservoir-modeling-concepts-and-equations) and [evaporation](https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/reservoir-modeling/reservoir-modeling-concepts-and-equations/evaporation). This project implements its own daily allocation/storage calculation, not HEC-HMS or Modified Puls routing. The reference supports mass conservation; the assumed operating rules and geometry are project scenarios.

The model is not calibrated or validated against an actual new structure. Siting, bathymetry, pumping energy, diversion permissions, subdaily floods, land availability, sedimentation, demands, channel losses and upstream operational feedback are not resolved. Positive numerical performance establishes conditional physical bookkeeping, not construction feasibility or an entitlement to divert.

## Reproduction and checks

```bash
python analysis/water_budget.py
python analysis/verify_water_budget.py
python -m unittest discover -s tests -v
```

No new packages beyond `requirements.txt` are required. The main run uses no random draws. Daily outputs are supplied for the baseline and nine nominal designs, with summary tables for the 324 combinations and 36 sensitivity checks. `annual_storage_ledgers.csv` includes 270 year-design budgets, and the summary stores input hashes. `verify_water_budget.py` checks saved budgets with Decimal arithmetic and separately reconstructs each day of the 10-million-m3, 10%-capture case directly from input CSVs. Unit checks cover empty storage, bypass protection, capacity limits, spill, transfer loss, carryover and mass conservation under varied inputs.
