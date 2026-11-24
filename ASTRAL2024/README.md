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

### Corrected-time MDV fallback

`recompute_mdv_fallback.jl` processes a date range with the corrected lidar
chunk index. The end date is exclusive. It requires the hourly stare NetCDFs,
`lidar_dt.jld2`, `file_beam_inds.jld2`, the sounding file under `data/netcdf/`,
and the PSL ship-navigation file. The examples below use its local path on
this machine; replace it if the file is stored elsewhere:

```sh
cd ASTRAL2024
micromamba run -n lidar julia --project=. recompute_mdv_fallback.jl \
  2024-05-19 2024-06-01 \
  /Users/deszoeks/Projects/ASTRAL/mast/data/SSTcompare/PSL/netcdf/ASTRAL-nav-met-sea-flux-10min_RV-Thompson_20240428_R1_thru_20240613.nc \
  epsilon_data/mdv_20240519_20240531.jld2
```

The script refuses to overwrite an existing output. Append a chunk ID to test
one chunk, for example `2516` with the date range `2024-05-20 2024-05-21`.
For an overlap comparison, append `--compare-vn`; aligned daily VectorNav
files must be present in `data/vn_sync_chunk_daily/`.

```sh
micromamba run -n lidar julia --project=. recompute_mdv_fallback.jl \
  2024-06-01 2024-06-10 \
  /Users/deszoeks/Projects/ASTRAL/mast/data/SSTcompare/PSL/netcdf/ASTRAL-nav-met-sea-flux-10min_RV-Thompson_20240428_R1_thru_20240613.nc \
  epsilon_data/mdv_vn_20240601_20240609.jld2 --compare-vn
```

The MDV is the mean over gates complete throughout each chunk; pitch and roll
are zero. Winds use the date-matched lidar product and interpolated soundings
minus ship speed, averaged equally when both exist, or the available source.
Lidar-wind and sounding levels are matched by gate index as in the production
workflow; their recorded height coordinates differ by 12 m. The fitted wind is
the chunk mean at each level. The shared TLS fit uses 48–125 m horizontal
separations and at least three fit bins. This range selection is not itself
proof that the flow is in the inertial range.

The output includes epsilon, binned structure functions, fit counts, MDV
vertical-velocity variance and sample counts, original chunk times and
GPS-clock-adjusted start times. Sentinel values are `-4`
for missing wind, `-3` for an insufficient fit and `-1` for negative slope.
With `--compare-vn`, `variance_mdv_paired`, `variance_vn_paired`,
`epsi_mdv_paired` and `epsi_vn_paired` use the same available beam samples,
wind and fit settings. `paired_samples` records the number of common samples.
Values remain missing when aligned motion is unavailable.

After processing the complete interval, render it instead of the provisional
shifted data with:

```sh
micromamba run -n lidar julia --project=. plot_tke_dissipation_may18_june09.jl \
  epsilon_data/mdv_20240519_20240531.jld2
```
