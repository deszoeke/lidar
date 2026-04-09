## Plan: QC Structure-Function Dissipation Fits

Add explicit fit-quality diagnostics to the structure-function regression used for A and dissipation, compute a simple composite confidence score (better than R^2 alone), and propagate diagnostics through chunk outputs so epsilon can be retained but tagged low-confidence.

**Steps**
1. Phase 1: Baseline and scope lock
2. Confirm the active computation path in /Users/deszoeks/Projects/lidar/ASTRAL2024/lidar_turbulence_cleanup.ipynb at the structure-function cell and in /Users/deszoeks/Projects/lidar/ASTRAL2024/lidar_turbulence_cleanup.jl for parity.
2.5 Copy the analysis to new files (for lidar_turbulence_qc-fit) to preserve the original.
3. Treat current fit model as fixed: D2 = noise + A*rho with epsilon derived from A.
4. Keep current physics inputs unchanged (pair construction, rho definition, equal-population binning, motion correction workflow).
5. Phase 2: Add per-level fit diagnostics where A is estimated
6. Extend regression logic in D2_rho_stare to compute and return all three requested diagnostics per level: R^2, normalized RMSE, and sample support (both raw pair count and binned-point count used in regression).
7. Add guard diagnostics for physically suspicious fits: non-positive slope A, non-finite outputs, negative noise, and unstable intercept/noise behavior.
8. Preserve existing A and noise outputs so downstream notebooks remain compatible.
9. Phase 3: Define confidence scoring (retain epsilon)
10. Implement a compact confidence score in [0,1] combining the three diagnostics (R^2, normalized RMSE, sample support), with simple monotonic transforms and clipping.
11. Keep score interpretable: 1 is excellent linear fit with adequate support; 0 is unusable fit.
12. There are 2 kinds of support, number of D2 pairs per bin, and number of bins with valid pairs used in the fit. Both should be saved and factored into confidence.
13. We don't a priori know thresholds for poor fits. Develop graphical review of epsilon vs confidence scores to identify appropriate cutoffs.
14. Phase 4: Propagate diagnostics through outputs
15. Update chunk-level arrays to store A, noise, R^2, normalized RMSE, support counts, and confidence score per height level alongside epsilon.
16. Persist these arrays in saved dissipation files so QC can be done post hoc without recomputing structure functions.
17. Keep current epsilon field and missing/sentinel conventions intact for backward compatibility.
18. Phase 5: QC visualization and review workflow
19. Add diagnostic plotting utilities for one chunk and campaign summaries: D2 vs rho with fitted line and residual behavior, vertical profiles of confidence metrics, and time-height maps of confidence.
20. Add quick selection helpers to compare epsilon distributions for all points vs high-confidence-only points.
21. Include simple failure counters (for example A<=0, low support, poor RMSE) to identify dominant failure modes.
22. Phase 6: Sensitivity checks
23. Evaluate threshold sensitivity for confidence classification (lenient vs strict) and report retained-data fraction vs fit quality.
24. Compare score behavior against raw R^2 to show why the composite score is more robust than variance-explained alone.

**Relevant files**
- /Users/deszoeks/Projects/lidar/ASTRAL2024/lidar_turbulence_cleanup.ipynb — structure-function cell and chunk processing sequence used interactively.
- /Users/deszoeks/Projects/lidar/ASTRAL2024/lidar_turbulence_cleanup.jl — D2_rho_stare, epsilon(A), and save pipeline used in script form.
- /Users/deszoeks/Projects/lidar/ASTRAL2024/read_epsilon.ipynb — downstream composite/QC analysis where confidence filters can be consumed.
- /Users/deszoeks/Projects/lidar/ASTRAL2024/epsilon_data — persisted per-chunk dissipation products to extend with QC fields.

**Verification**
1. Unit-style check on a single chunk: verify returned arrays have expected shapes and finite values where bins are valid.
2. Regression sanity check at selected heights: fitted D2 values track binned D2 with expected residual magnitude, and metrics change appropriately when fit degrades.
3. Confidence behavior check: deliberately poor-support or noisy levels produce lower confidence than clean levels.
4. Persistence check: save and reload one output file and confirm confidence/diagnostic fields round-trip with epsilon.
5. Downstream check: filter by confidence in read_epsilon workflow and confirm resulting composites change in expected direction (reduced tails, improved fit consistency).

**Decisions**
- Include all three metrics: R^2, normalized RMSE, and sample support.
- Keep epsilon values even when fit is poor; encode confidence separately.
- Use a composite confidence score rather than R^2 alone to avoid over-trusting sparse or noisy fits.
- Avoid changing core structure-function physics at this stage; focus on QC observability and scoring.

**Further Considerations**
1. Score design recommendation: weighted harmonic-style combination of normalized R^2, inverse RMSE term, and support term to strongly penalize any single weak dimension.
2. Backward compatibility recommendation: add new fields without renaming existing epsilon keys so old notebooks continue to run.
3. Optional later extension: introduce confidence-weighted composites in mixed-layer analysis instead of hard confidence thresholds.