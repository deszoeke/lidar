# Production Run Additions

This note records the main production-run work added in the April 2026 notebook updates to `lidar_turbulence_production.ipynb`.

## Principles of calculation

The production retrieval is based on a linear fit of the binned structure function at each height level,

$$
D^2(\rho) = A\rho + \mathrm{noise}.
$$

The fitted slope `A` is converted to turbulent kinetic energy dissipation rate using the existing longitudinal structure-function convention,

$$
\epsilon(A) = \left(\sqrt{\tfrac{3}{4}A/C_{2ll}}\right)^3,
$$

with `C2ll = 2.0`.

The fit-quality additions were kept non-gating. The production epsilon path remains:
- `A` missing or too few usable bins: `epsi = -3`
- `A < 0`: `epsi = -1`
- missing wind for the chunk: `epsi = -4`
- not yet computed: `epsi = -5`
- otherwise: `epsi = epsilon(A)`

The uncertainty work stores slope uncertainty directly and leaves epsilon uncertainty as a postproduction calculation. The first-order propagation used later is

$$
\mathrm{se}_\epsilon = \tfrac{3}{2}\,\epsilon(A)\,\frac{\mathrm{se}_A}{A},
$$

for valid `A > 0` and finite `se_A`.

## Efficiency changes

The main efficiency change was replacing repeated fit-stat passes with a single helper, `fit_stats_onepass(x, y)`, applied once per height level after binning.

That helper:
- extracts valid finite `x, y` only once
- centers the vectors once
- reuses the same sums to compute slope, residual variance, `R2`, `RMSE`, `se_A`, and confidence intervals

In practice this avoids calling the valid-data extraction separately for each metric. The notebook discussion during this work was explicitly to avoid recomputing the same masked fit inputs three or more times per level.

The implementation also makes `nbins` exact: it is the number of binned regression points actually used in the fit, not a raw pair count.

## Output variables saved

The production notebook now saves the following variables to the JLD output:
- `epsi`
- `dtime_st`
- `dtime_en`
- `height`
- `A`
- `se_A`
- `nbins`
- `R2_tmp`
- `RMSE_tmp`
- `A_lo_tmp`
- `A_hi_tmp`
- `epsi_lo_tmp`
- `epsi_hi_tmp`

The intent of these fields is:
- `epsi`: production dissipation retrieval with sentinel coding
- `A`: fitted structure-function slope at each chunk and height
- `se_A`: standard error of the fitted slope
- `nbins`: number of binned fit points used in the regression
- `R2_tmp`, `RMSE_tmp`: non-gating fit diagnostics
- `A_lo_tmp`, `A_hi_tmp`: 95% confidence interval bounds on slope
- `epsi_lo_tmp`, `epsi_hi_tmp`: epsilon interval bounds obtained by transforming the slope CI through `epsilon(A)`

## Postproduction utility

The notebook also now includes a tiny loader utility that reads a saved JLD file, loads `A` and `se_A`, and computes:
- `se_epsilon`
- `rel_se_epsilon = se_epsilon / epsi`

This loader accepts both canonical saved names (`A`, `se_A`) and fallback temporary names (`A_tmp`, `se_A_tmp`) for compatibility with older outputs.