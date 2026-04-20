# Implementation Plan: Fix Double-Counting Bug & Add Prior Reset

## Overview
This plan fixes the double-counting of `prior_offset` in the offset calculation pipeline and adds configurable prior reset modes to prevent cascade failures from bad syncs.

---

## Bug Fix: Double-Counting of `prior_offset`

### Root Cause
- **Line 764** in `refine_offset_20hz()` computes: `final_total = prior_offset + final_native`
- **Line 869** in `process_sync_data()` re-adds prior: `final_native = ... r.prior_offset + final_native_resid`
- Result: `prior_offset` is added twice, producing incorrect stored totals

### Current Behavior (BROKEN)
```
prior_offset = -1.02 s
xcorr20.lag_seconds = -1.0 s

refine_offset_20hz returns:
  final_offset_native = -1.02 + (-1.0) = -2.02 s  ← Already includes prior!

process_sync_data stores:
  final_native = prior_offset + final_native_resid 
               = -1.02 + (-2.02)  ← BUG! Adds prior again
               = -3.04 s  ← Wrong!
```

### Correct Behavior (FIXED)
```
prior_offset = -1.02 s
xcorr20.lag_seconds = -1.0 s

refine_offset_20hz returns:
  final_offset_native = -1.0 s  ← Residual only, NOT prior + residual

process_sync_data stores:
  final_native = prior_offset + final_native_resid 
               = -1.02 + (-1.0)
               = -2.02 s  ← Correct!
```

---

## Changes Required

### Change 1: Fix `refine_offset_20hz()` to return residual only

**File**: `lidar_vn_sync.jl`, Lines ~764–766

**Before:**
```julia
    final_native = round(xcorr20.lag_seconds / dt20) * dt20
    final_round_s = round(Int, final_native)
    final_total = prior_offset + final_native
```

**After:**
```julia
    final_native = round(xcorr20.lag_seconds / dt20) * dt20
    final_round_s = round(Int, final_native)
    # NOTE: final_native is the 20 Hz residual ONLY (not prior + residual)
    # Caller (process_sync_data) will compose with prior_offset separately
```

**Reasoning**: 
- `final_native` should be the 20 Hz residual relative to the 1 Hz stage
- The caller `process_sync_data()` handles prior composition
- This makes offsets composable at each stage: `total = prior + residual_1hz + residual_20hz`

### Change 2: Update `process_sync_data()` to not re-add prior

**File**: `lidar_vn_sync.jl`, Lines ~869–870

**Before:**
```julia
    final_native_resid = offset_or_sentinel(rf.final_offset_native; sentinel=offset_sentinel_seconds)
    final_native = isfinite(final_native_resid) && isfinite(r.prior_offset) ? r.prior_offset + final_native_resid : offset_sentinel_seconds
```

**After:**
```julia
    final_native_resid = offset_or_sentinel(rf.final_offset_native; sentinel=offset_sentinel_seconds)
    # final_offset_native is already the 20 Hz residual; just add to 1 Hz result (which includes prior)
    final_native = isfinite(final_native_resid) && isfinite(s1.final_offset) ? s1.final_offset + final_native_resid : offset_sentinel_seconds
```

**Reasoning**:
- `s1.final_offset` already includes `prior_offset` (it's `prior + coarse + fine`)
- `final_native_resid` is now the 20 Hz residual (no prior in it)
- Composition is correct: `s1.final_offset + residual_20hz`

### Change 3: Update fallback result in `process_sync_data()` to not add prior

**File**: `lidar_vn_sync.jl`, Lines ~851–853

**Before:**
```julia
    push!(seq_ref20, (; r..., lidar_dt=w.stare_dt, final_offset_20hz=offset_sentinel_seconds, final_offset_native=offset_sentinel_seconds, final_offset_round_s=round(Int, offset_sentinel_seconds), vn2_1s_aligned=fill(NaN, length(w.stare_dt)), pitch_1s_aligned=fill(NaN, length(w.stare_dt)), roll_1s_aligned=fill(NaN, length(w.stare_dt)), mdv_residual_1s=fill(NaN, length(w.stare_dt))))
```

**Note**: In the case where `!isfinite(rf.final_offset_20hz)`, the fallback should be to use 1 Hz result directly:

**Update** (surrounding context, lines ~845–853):
```julia
else  # rf.final_offset_20hz is NaN, use 1 Hz result directly
    final_20 = offset_sentinel_seconds
    final_native_resid = offset_sentinel_seconds
    final_native = isfinite(s1.final_offset) ? s1.final_offset : offset_sentinel_seconds
    final_round = isfinite(final_native) ? round(Int, final_native) : round(Int, offset_sentinel_seconds)
    
    push!(seq_ref20, (; r..., lidar_dt=w.stare_dt, final_offset_20hz=final_20, final_offset_native=final_native, final_offset_round_s=final_round, vn2_1s_aligned=fill(NaN, length(w.stare_dt)), pitch_1s_aligned=fill(NaN, length(w.stare_dt)), roll_1s_aligned=fill(NaN, length(w.stare_dt)), mdv_residual_1s=fill(NaN, length(w.stare_dt))))
end
```

---

## Enhancement: Configurable Prior Reset

### Motivation
Currently, if one chunk's sync fails or produces a wrong offset, that bad offset carries forward as prior to all subsequent chunks in the file, poisoning the entire file's estimates. 

**Solution**: Add configurable reset modes so users can:
1. Carry prior forward (default, current behavior)
2. Reset prior at file boundaries
3. Reset prior when a sync is rejected (low confidence)

### Change 4: Add reset parameters to `process_sync_data()`

**File**: `lidar_vn_sync.jl`, Lines ~776–785 (function signature)

**Before:**
```julia
function process_sync_data(beams, env, Vn, UV, ic_list; 
    ntop=80, jump_threshold_seconds=0.8, min_fine_peak=0.08, 
    backward_windows=2, backward_tol=0.35, max_gap_samples=3,
    min_accept_peak_norm=0.08, fallback_search_seconds=60.0,
    offset_sentinel_seconds=-9999.0, offset_sentinel_ms=Int64(-9_999_000),
    nan_log_path=joinpath("epsilon_data", "nan_offset_chunks.log"), 
    reset_nan_log=true, 
    iter_log_path=joinpath("epsilon_data", "vn_log_$(Dates.format(Dates.now(), dateformat"yyyymmdd_HHMMSS")).txt"), 
    reset_iter_log=true )
```

**After:**
```julia
function process_sync_data(beams, env, Vn, UV, ic_list; 
    ntop=80, jump_threshold_seconds=0.8, min_fine_peak=0.08, 
    backward_windows=2, backward_tol=0.35, max_gap_samples=3,
    min_accept_peak_norm=0.08, fallback_search_seconds=60.0,
    offset_sentinel_seconds=-9999.0, offset_sentinel_ms=Int64(-9_999_000),
    nan_log_path=joinpath("epsilon_data", "nan_offset_chunks.log"), 
    reset_nan_log=true, 
    iter_log_path=joinpath("epsilon_data", "vn_log_$(Dates.format(Dates.now(), dateformat"yyyymmdd_HHMMSS")).txt"), 
    reset_iter_log=true,
    reset_prior_on_file_boundary::Bool=false,  # NEW
    reset_prior_on_rejected_sync::Bool=false   # NEW
)
```

### Change 5: Pass reset parameters to `run_sequential_offsets()`

**File**: `lidar_vn_sync.jl`, Lines ~816–817

**Before:**
```julia
    seq_results = run_sequential_offsets(beams, env, Vn, UV, ic_list; 
        ntop=ntop, jump_threshold_seconds=jump_threshold_seconds, ...
```

**After:**
```julia
    seq_results = run_sequential_offsets(beams, env, Vn, UV, ic_list; 
        ntop=ntop, jump_threshold_seconds=jump_threshold_seconds, 
        ...,
        reset_prior_on_file_boundary=reset_prior_on_file_boundary,
        reset_prior_on_rejected_sync=reset_prior_on_rejected_sync
```

### Change 6: Add reset logic in `run_sequential_offsets()`

**File**: `lidar_vn_sync.jl`, Lines ~450–550 (function signature and body)

Add to function signature:
```julia
function run_sequential_offsets(beams, env, Vn, UV, ic_list; 
    ...existing params...,
    reset_prior_on_file_boundary::Bool=false,
    reset_prior_on_rejected_sync::Bool=false
)
```

Add file tracking and reset logic in main loop (around line 470–490):
```julia
    results = NamedTuple[]
    ifile_prev = -1
    
    for ic in ic_list
        # Determine current file
        ifile = findfirst(ic .<= env.bigind_file_ien)
        
        # Extract sync window
        w = extract_sync_window(beams, env, state_ref, Vn, UV, ic; ntop=ntop, nc_dir=nothing)
        
        # Compute prior, with optional resets
        if reset_prior_on_file_boundary && ifile != ifile_prev
            # Reset prior at file boundary
            prior_from_result = prior_from_history(results, w.stare_dt[1]; reset_to=0.0)
        elseif reset_prior_on_rejected_sync && !isempty(results) && !get(results[end], :accepted_sync, true)
            # Reset prior if last sync was rejected
            prior_from_result = prior_from_history(results, w.stare_dt[1]; reset_to=0.0)
        else
            # Normal carry-forward
            prior_from_result = prior_from_history(results, w.stare_dt[1])
        end
        
        ifile_prev = ifile
        
        # ... rest of loop continues as before, using prior_from_result ...
```

### Change 7: Update `prior_from_history()` signature to support forced resets

**File**: `lidar_vn_sync.jl`, Lines ~550–560 (function signature)

**Before:**
```julia
function prior_from_history(results, start_dt)
```

**After:**
```julia
function prior_from_history(results, start_dt; reset_to::Union{Float64, Nothing}=nothing)
```

**In function body** (early lines):
```julia
    if !isnothing(reset_to)
        return (; prior_seconds=reset_to, previous_offset=NaN, recent_offset=NaN, hourly_offset=NaN, weight_prev=0.0, weight_recent=0.0, weight_hourly=0.0)
    end
    
    # ... rest of normal prior_from_history logic ...
```

---

## Testing Plan

### Unit Test 1: Verify double-count fix
**Location**: Notebook cell or separate test file

**Test code**:
```julia
# Simulate the offset composition
prior_s = -1.02
residual_20hz = -1.0

# Old way (WRONG):
old_total = prior_s + (prior_s + residual_20hz)  # = -2.02 - 1.0 = -3.04
@test old_total == -3.04  # Confirm bug

# New way (FIXED):
new_total = prior_s + residual_20hz  # = -1.02 + (-1.0) = -2.02
@test new_total ≈ -2.02  # Should match stored_total from accounting trace
```

### Integration Test 2: Re-run accounting trace after fix
**Location**: Notebook cell (add after fix implementation)

**Test code** (add to notebook):
```julia
# Re-run accounting trace with fixed code
# Expected: stored_total should match best_grid_shift or be closer to it
# Old trace showed: stored_total = -2.02, best_grid_shift = -0.56 to -0.9
# New trace should show: stored_total ≈ best_grid_shift
```

### Integration Test 3: Test prior reset modes
**Location**: Notebook cell with synthetic test chunks

**Test code**:
```julia
# Create test chunks spanning a file boundary
# Run with reset_prior_on_file_boundary=true
# Verify: prior resets to 0.0 when file changes
# Verify: prior carries forward when false (default)
```

### Regression Test 4: Compare old vs new results
**Location**: Full pipeline test

**Test code**:
```julia
# Run process_sync_data on full dataset with:
# 1. reset_prior_on_file_boundary=false (default, should match old behavior mostly)
# 2. reset_prior_on_file_boundary=true (new safe behavior)
# Verify: flag=true produces different but sensible results
```

---

## Deployment Steps

1. **Make Changes 1–3** (bug fix)
   - Modify `refine_offset_20hz()` line 764
   - Modify `process_sync_data()` line 869
   - Modify fallback at line 851

2. **Test Changes 1–3**
   - Run Unit Test 1 (confirm -2.02 composition)
   - Run Integration Test 2 (re-run accounting trace)
   - Verify `stored_total` no longer double-counts

3. **Make Changes 4–7** (prior reset enhancement, optional)
   - Add parameters to `process_sync_data()` signature
   - Add reset logic to `run_sequential_offsets()`
   - Update `prior_from_history()` signature

4. **Test Changes 4–7**
   - Run Integration Test 3 (file boundary reset)
   - Run Regression Test 4 (full pipeline with both modes)

5. **Update Documentation**
   - Add docstring to new parameters
   - Reference OFFSET_EQUATIONS_AND_RESET_STRATEGY.md in module comments

---

## Backward Compatibility

- **Changes 1–3** are fixes that correct the math; they're not backwards-compatible but necessary
- **Changes 4–7** are enhancements with default `false`, so old code continues to work
- All changes are isolated to `lidar_vn_sync.jl` (no changes needed to readers, callers, or saved data)

---

## Estimated Impact

- **Lines changed**: ~10–15 lines in `lidar_vn_sync.jl`
- **Functions affected**: `refine_offset_20hz()`, `process_sync_data()`, `run_sequential_offsets()`, `prior_from_history()`
- **Test coverage needed**: 4 integration tests (covered in Testing Plan)
- **Risk level**: Low (bug fix is isolated, enhancements are opt-in)
