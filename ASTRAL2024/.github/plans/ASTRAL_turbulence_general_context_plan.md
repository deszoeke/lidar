# ASTRAL Turbulence: General Context + Plan

## Context (What Is True)
- Campaign: ASTRAL in the Bay of Bengal during monsoon onset.
- Core idea: infer boundary-layer turbulence structure by combining shipboard lidar dissipation profiles with surface flux forcing.
- Key variables:
  - epsilon(z,t): turbulent kinetic energy dissipation rate from lidar-derived products
  - B0(t): surface buoyancy flux from PSL met/flux data
  - h_mix(t): mixing-layer height inferred from vertical epsilon structure
  - Monin-Obukhov length:
  $L = - {u^*}^3 / (k B_0)$
  (by convention <0 for unstable conditions $B_0>0$), the vertical scale at which the buoyancy integral equals the mechanical turbulence generation.
- Typical analysis challenge: distributions are skewed/intermittent, so median and mean can imply different physics.

## Plan (What To Do Next)
1. Ingest and align datasets in time:
   - epsilon profiles, profile timestamps, PSL flux variables.
2. Compute and QC forcing terms:
   - derive B0, mask missing/non-finite values, enforce time-match tolerance.
3. Compute structure metrics:
   - estimate h_mix, then derive normalized coordinates (e.g., z/h_mix).
4. Build composites:
   - conditional by regime (e.g., sign/magnitude of B0, wind bins, rain flags).
5. Report robust statistics:
   - sample count, median, mean, IQR, upper-tail quantiles.
6. Interpret physically:
   - separate background buoyancy-driven mixing from intermittent mechanical/shear-driven events.

## Rule Of Thumb: Context vs Plan
- Context = assumptions, definitions, data provenance, and known constraints.
- Plan = ordered actions you intend to execute next.
- Keep both: context prevents drift; plan drives progress.
