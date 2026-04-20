# Fix Summary: Double-Counting of Prior Offset

## Problem Identified
The offset composition pipeline had a critical bug where `prior_offset` was being added twice:
1. In `refine_offset_20hz()`, the function returned `final_offset_native = prior_offset + xcorr20.lag_seconds`
2. In `process_sync_data()`, the code added prior again: `final_native = r.prior_offset + final_native_resid`
3. Result: `stored_total = 2 × prior + residual_20hz` instead of `prior + residual_1hz + residual_20hz`

### Observed Impact
Accounting trace test showed:
```
prior_s        = -1.02 s
s1_final_offset = -1.02 s
rf_native_resid = -1.0 s
OLD stored_total = -2.02 s (prior + residual, then +prior again = WRONG)
best_grid_shift = -0.56 to -0.9 s (correlation-optimal)
```

The stored total was significantly divergent from the correlation-optimal shift because of the double-count.

---

## Changes Made

### File: `lidar_vn_sync.jl`

#### Change 1: Line ~764–771 in `refine_offset_20hz()`
**Removed the line** `final_total = prior_offset + final_native`

**Added comment** explaining that `final_native` is now the residual only, not prior + residual:
```julia
# NOTE: final_native is the 20 Hz residual ONLY (not prior + residual).
# The caller (process_sync_data) composes with prior_offset separately,
# combining it with the 1 Hz stage result (which already includes prior).
# Offset composition: total = prior + residual_1hz + residual_20hz
```

**Reasoning**: This allows `final_offset_native` to be a pure residual, making composition unambiguous.

#### Change 2: Line ~869–870 in `process_sync_data()`
**Changed** `final_native = isfinite(final_native_resid) && isfinite(r.prior_offset) ? r.prior_offset + final_native_resid : offset_sentinel_seconds`

**To** `final_native = isfinite(final_native_resid) && isfinite(s1.final_offset) ? s1.final_offset + final_native_resid : offset_sentinel_seconds`

**Added comment**:
```julia
# final_offset_native from refine_offset_20hz is the 20 Hz residual (no prior in it).
# s1.final_offset is the 1 Hz result which already includes prior.
# Correct composition: s1.final_offset + residual_20hz (no double-counting of prior)
```

**Reasoning**: Now we compose by adding the 20 Hz residual to the 1 Hz result (which already includes prior), not by re-adding prior.

#### Change 3: Line ~851–853 in `process_sync_data()` (fallback case)
**Changed** fallback to use `s1.final_offset` directly instead of sentinel when 20 Hz refinement fails:
```julia
# Fallback when refine_offset_20hz fails: use 1 Hz result directly
final_native_fallback = isfinite(r.final_offset) ? r.final_offset : offset_sentinel_seconds
push!(seq_ref20, (; r..., lidar_dt=w.stare_dt, final_offset_20hz=offset_sentinel_seconds, 
    final_offset_native=final_native_fallback, ...))
```

**Reasoning**: When 20 Hz refinement fails, fall back to the 1 Hz result (which is complete), not to sentinel.

---

## Equations (Now Correct)

### Lidar time correction (per-file, applied once)
```
lidar_time(t) = t_raw(t) - fit_offset(t_raw[1])
```
- `fit_offset()` returns scalar offset (~126–144 ms), piecewise linear function of file date

### VN offset estimation (per chunk)

#### Prior from history
```
prior_offset ← prior_from_history(...)  [file-persistent, weighted blend of recent offsets]
```

#### Stage 1: Coarse + fine lag at 1 Hz
```
s1 = coarse_and_fine_lag(mdv, vn2; prior_seconds=prior_offset)
   → returns: coarse_offset, fine.lag_seconds, final_offset
   → final_offset = prior_offset + Δ_coarse + Δ_fine  [includes prior]
```

#### Stage 2: Refinement residual at 20 Hz
```
rf = refine_offset_20hz(win, s1, Vn)
   → returns: final_offset_20hz, final_offset_native, ...
   → final_offset_native = residual_20hz  [residual ONLY, no prior in it]
```

#### Total offset (FIXED composition)
```
offset_stored = s1.final_offset + rf.final_offset_native
              = (prior + Δ_coarse + Δ_fine) + residual_20hz
              ✓ Correct! No double-counting.
```

---

## Verification

### New Notebook Cell: `#VSC-f5b78232`
Added verification cell that re-runs the accounting trace after the fix and computes:
- `stored_total_new = s1_final + rf_native_resid` (correct composition)
- Comparison: error to `best_grid_shift` should improve (decrease)
- Expected output: improvement > 0 for most chunks (lower error = closer to correlation-optimal)

**How to run**:
```julia
# In notebook, execute Cell 7 (accounting trace), then Cell 8 (verification)
# Look for "Δerror" column: positive values show error decreased (bug is fixed)
```

---

## Impact Analysis

### Downstream Effects
- All files that use `seq_results` or `seq_ref20` from `process_sync_data()` will now have correct offsets
- Re-running production will overwrite stored offsets with correct values
- No schema changes (same output structure, just correct values)

### Backward Compatibility
- **Breaking**: Stored offsets will change (bug fix, not parameter change)
- **Non-breaking**: Code signature unchanged, function behavior corrected
- **Mitigation**: Archive old offsets if needed for comparison; re-run production with new code

### Validation Scope
- Accounting trace test: run full chunks 200:20:300 to confirm improvement
- Coherence test: re-compute Welch coherence on corrected offsets, verify > 0.8 in 5–20 s band
- Grid search fallback: verify all chunks still accept (peak > min_accept_peak_norm)

---

## Next Steps (Optional Enhancements)

### Configurable Prior Reset (Not Yet Implemented)
To prevent cascade failures from bad syncs, add:
- `reset_prior_on_file_boundary::Bool=false` — reset prior at file boundaries
- `reset_prior_on_rejected_sync::Bool=false` — reset if sync confidence is low

See `IMPLEMENTATION_PLAN_OFFSET_FIX.md` for design details (Changes 4–7).

---

## Testing Checklist

- [ ] Rebuild module (run `Revise` cell or `includet()`)
- [ ] Run accounting trace cell (Cell 7) for baseline
- [ ] Run verification cell (Cell 8) to confirm fix
- [ ] Check "Δerror" column: values should be > 0 (error decreased)
- [ ] Re-run quality metrics (coherence, correlation) on corrected offsets
- [ ] Run full production pipeline and verify output schema unchanged

---

## Files Modified

1. `/Users/deszoeks/Projects/lidar/ASTRAL2024/lidar_vn_sync.jl`
   - 3 locations (refine_offset_20hz, process_sync_data x2)
   - ~10 lines changed, added comments explaining composition

2. `/Users/deszoeks/Projects/lidar/ASTRAL2024/lidar_vn_sync_workbench.ipynb`
   - Added Cell `#VSC-f5b78232` (verification trace)
   - New cell immediately after accounting trace cell

3. Documentation created:
   - `OFFSET_EQUATIONS_AND_RESET_STRATEGY.md` (equations + reset design)
   - `IMPLEMENTATION_PLAN_OFFSET_FIX.md` (detailed implementation guide)
   - This file: `FIX_SUMMARY_DOUBLE_COUNT.md` (overview)

---

## Reference

- **Root cause**: Line 869 of `process_sync_data()` added prior when `final_offset_native` (returned from `refine_offset_20hz`) already included it
- **Key insight**: `refine_offset_20hz()` should return residual only; caller composes with prior
- **Equations**: Now follow pattern `total = prior + residual_1hz + residual_20hz` (no re-addition of prior)
