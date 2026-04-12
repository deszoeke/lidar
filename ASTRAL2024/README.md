# Routines for reading and processing Doppler wind lidar 
for the 2024 ASTRAL/EKAMSAT cruise.

(c) Simon de Szoeke

## Data to read
  - Halo Photonics XR:
    `User` wind profiles and
    `Stare` vertical stares
  - VectorNav accelerometer
  - Navigation: ship's navigation message or NOAA PSL  

## Turbulence processing
Estimate vertical velocity variance and tubulent kinetic energy (TKE) dissipation rate from vertical stares.

### Production workflow
- The main production workflow is implemented in `lidar_turbulence_production.ipynb`.
- Inputs are Halo StreamLineXR vertical stare beams, collocated VectorNav motion, and mean relative wind from the EKAMSAT netCDF product.
- Data are processed chunk-by-chunk over the stare record, loading lidar files on demand into a periodic buffer and restricting production to chunks with overlapping VectorNav coverage.
- For each height level in a chunk, the code forms a structure function in binned separation space and fits
  $D^2(\rho) = A\rho + \mathrm{noise}$
  where $A$ is the slope used to retrieve dissipation.
- Dissipation is derived from the fitted slope using the existing longitudinal structure-function convention
  $\epsilon(A) = \left(\sqrt{\tfrac{3}{4}A/C_{2ll}}\right)^3$
  with `C2ll = 2.0`.
- Production keeps the epsilon retrieval behavior unchanged and stores fit-quality quantities as diagnostic outputs only; no additional QC gating is applied during retrieval.

### Saved production outputs
- The full production run writes a JLD file in `epsilon_data/` named `epsi_stare_chunks_rYYYYMMDD_HHMM.jld`.
- Persisted variables are `epsi`, `dtime_st`, `dtime_en`, `height`, `A`, `se_A`, `nbins`, `R2_tmp`, `RMSE_tmp`, `A_lo_tmp`, `A_hi_tmp`, `epsi_lo_tmp`, and `epsi_hi_tmp`.
- `A`, `se_A`, and `nbins` are the canonical saved fit outputs. Here `nbins` means the number of binned regression points used in the fit, not the number of raw sample pairs.
- `epsi` uses negative sentinel values to retain production-state information: `-5` uncomputed default, `-4` missing wind, `-3` missing or insufficient fit, and `-1` negative fitted slope.

### Fit statistics and postprocessing
- Recent production updates compute fit diagnostics in one pass per height level from the same valid binned `rho` and `D2` samples used for the slope fit.
- The one-pass helper reuses centered sums to derive `A`, `R2`, `RMSE`, `se_A`, and 95% confidence intervals without repeatedly re-filtering the same binned data.
- Confidence intervals on epsilon are obtained by transforming the slope confidence interval endpoints through `epsilon(A)`.
- A postproduction loader utility in the notebook computes `se_epsilon` and relative uncertainty `se_epsilon / epsi` from saved `A` and `se_A`.
