# Timing Offset Equations & Prior Reset Strategy

## Current Implementation (Broken)

### Equation 1: Lidar time synchronization
```
lidar_time(t) = t_lidar_raw(t) - lidar_clock_fast_by(t_file_start)
```
- `lidar_clock_fast_by` = `fit_offset(stare_dt_raw[1])` per file (scalar, ~126–144 ms)
- Applied once per file in `chunk_lidar_datetimes()`

### Equation 2: 1 Hz coarse + fine lag stage
```
s1_result = coarse_and_fine_lag(mdv, vn2; prior_seconds)
  ├─ coarse_offset = prior_seconds + iterative_coarse_search(...)
  └─ fine.lag_seconds = fft_xcorr_lag(...) relative to coarse_offset
  
final_offset_1hz = s1_result.final_offset 
                 = s1_result.coarse_offset + s1_result.fine.lag_seconds
                 = prior_seconds + Δ_coarse + Δ_fine
```

### Equation 3: 20 Hz refinement stage
```
rf_result = refine_offset_20hz(win, s1_result, Vn)
  → uses prior_offset = s1_result.prior_seconds to set VN window

rf_result.final_offset_native = prior_offset + xcorr20.lag_seconds
```

### Equation 4: What gets stored (BUG)
```
stored_total = prior_offset + final_offset_native
             = prior_offset + (prior_offset + xcorr20.lag_seconds)
             = 2 × prior_offset + xcorr20.lag_seconds  ← DOUBLE COUNT!
```

## Observed Behavior
Trace (chunks 220–280):
| Field | Value | Source |
|-------|-------|--------|
| `prior_s` | `-1.02` | from prior history |
| `s1_final_offset` | `-1.02` | coarse centered on prior + fine, got -1.02 again (stuck at prior) |
| `rf_native_resid` | `-1.0` | 20 Hz refinement residual relative to prior |
| `stored_total` | `-2.02` | prior + residual ← **double counts prior** |
| `best_grid_shift` | `-0.56 to -0.9` | correlation-optimal, not returned |

**Root cause**: `refine_offset_20hz` returns `final_offset_native = prior_offset + xcorr_lag`, but then `process_sync_data` computes `final_native = prior_offset + final_offset_native`, creating a double-count.

---

## Correct Equations (What We Want)

### Stage 1: Coarse + fine at 1 Hz, residual relative to prior
```
Δ_coarse_fine = coarse_and_fine_lag(mdv, vn2; prior_seconds=prior_offset)
              → returns offset relative to prior_offset
```

### Stage 2: Refinement residual at 20 Hz
```
rf = refine_offset_20hz(win, Δ_coarse_fine, Vn)
   → returns final_offset_native = residual relative to Δ_coarse_fine
   
total_offset_chunk = prior_offset + Δ_coarse_fine + final_offset_20hz_residual
```

### Stage 3: Carry forward or reset
```
prior_offset_next = (accept_sync) ? total_offset_chunk : prior_offset
                  OR
                    (reset_at_file_boundary) ? 0.0 : prior_offset
                  OR
                    (reset_on_failed_sync) ? 0.0 : prior_offset
```

---

## Proposed API for Prior Reset (Minimal, Safe)

Add to `process_sync_data()` signature:
```julia
function process_sync_data(
    beams, env, Vn, UV, ic_list; 
    ...existing params...,
    prior_reset_mode::Symbol=:carry_forward,  # NEW
    reset_on_file_boundary::Bool=false,       # NEW
    reset_on_failed_sync::Bool=false,         # NEW
)
```

### Reset modes:
1. **`:carry_forward`** (default) – inherit prior from previous chunk, even across files
2. **`:reset_each_file`** – set `prior = 0.0` at file boundaries
3. **`:reset_on_bad_quality`** – set `prior = 0.0` if `accepted_sync == false`

### Implementation location (4 changes in lidar_vn_sync.jl):

#### Change 1: `run_sequential_offsets()` signature
Add parameter `reset_on_failed_sync::Bool=false` and pass to `prior_from_history()`.

#### Change 2: `prior_from_history()` logic
```julia
function prior_from_history(history, start_dt, ifile, ifile_prev; 
    reset_on_file_boundary::Bool=false, reset_on_failed_sync::Bool=false)
    
    # Reset at file boundary
    if reset_on_file_boundary && ifile != ifile_prev
        return (; prior_seconds=0.0, previous_offset=NaN, ...)
    end
    
    # Reset if last sync was rejected
    if reset_on_failed_sync && !isempty(history) && !get(history[end], :accepted_sync, true)
        return (; prior_seconds=0.0, previous_offset=NaN, ...)
    end
    
    # Normal carry-forward logic
    ...
end
```

#### Change 3: Track file index in `run_sequential_offsets()`
```julia
ifile_current = 0
for (idx, ic) in enumerate(ic_list)
    ifile_needed = findlast(env.bigind_file_starts .<= env.ists[ic])
    prior = prior_from_history(history, start_dt, ifile_needed, ifile_current; 
        reset_on_file_boundary, reset_on_failed_sync)
    ifile_current = ifile_needed
    ...
end
```

#### Change 4: Fix the double-count bug in `process_sync_data()`
```julia
# Current (WRONG):
final_native = r.prior_offset + rf.final_offset_native

# Fixed:
final_native = s1.final_offset + rf.final_offset_native
             # (s1.final_offset already includes prior_offset)
```

---

## Validation Plan

1. **Unit test**: Verify `prior_from_history()` returns 0.0 when file boundary crossed with `reset_on_file_boundary=true`
2. **Integration test**: Run `process_sync_data()` on chunks spanning a file boundary, verify prior does/doesn't carry per mode
3. **Regression test**: Run with `prior_reset_mode=:carry_forward` (default), verify prior carries as before
4. **Smoke test**: Re-run accounting trace on chunks 220–280 with each mode, verify:
   - `stored_total` no longer double-counts
   - `stored_total` approaches `best_grid_shift` when using fallback-search or wide coarse window

---

## Safety Notes

- **Backward compatible**: default `reset_on_file_boundary=false` preserves current (buggy) behavior for now
- **Isolated to sync module**: file tracking stays in `run_sequential_offsets()`, no reader changes needed
- **Gradual migration**: user can opt in per call without rebuilding pipelines
