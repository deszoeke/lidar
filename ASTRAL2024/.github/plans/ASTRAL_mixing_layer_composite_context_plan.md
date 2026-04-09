# ASTRAL Mixing-Layer Composite: Session Context + Plan

## Context (Current Analysis State)
- Objective: composite epsilon normalized by surface forcing on stretched coordinate z/h_mix.
- Current composite setup:
  - eta in [0,2] with 40 bins (20 below h_mix, 20 above).
  - epsilon filter: keep finite, positive epsilon values.
  - h_mix from epsilon profile crossing criterion (computed separately).
  - B0 matched to profile time using nearest valid PSL timestamp with guard window.
- Regime split:
  - requested separate composites for B0>0 and B0<0, without abs(B0) in normalization.
  - in current matched period, B0<0 appears absent (or negligible after QC).
- Outlier handling:
  - weak-flux filter applied: exclude profiles with |B0| < 1e-4.
  - distribution remains right-skewed; median << mean due to intermittent high values.

## Plan (Immediate Next Steps)
1. Lock a baseline figure set:
   - unstable composite (B0>0) with median, mean, IQR, and counts.
2. Add robustness checks:
   - sensitivity to |B0| threshold (1e-4, 2e-4, 5e-4).
   - compare linear vs log-x display of epsilon/B0.
3. Quantify intermittency explicitly:
   - compute P90/P95/P99 by eta bin and contribution of top decile to the mean.
4. Add physical conditioning:
   - stratify by wind speed and precipitation flag to isolate mechanical bursts.
5. Produce interpretation-ready summary:
   - concise statement of typical state (median) vs event-driven mean enhancement.

## Minimal Decision Log
- Use median + IQR for typical structure.
- Use mean + high quantiles to represent intermittent energetic events.
- Keep context and plan separate in notes to reduce confusion and prompt drift.
