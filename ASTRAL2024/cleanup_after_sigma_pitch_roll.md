# Cleanup Plan After Pitch/Roll Uncertainty Implementation

## Test Organization

### Issue: Tests clutter notebook flow
**Current:** TEST 1-3 cells in lidar_turbulence_production.ipynb (cells 23-25)

**Proposal:**
- Move TEST 2-3 (displacement_variance_pitch_roll, propagate_rho_uncertainty tests) to a testing module in `displacement_uncertainty.jl`
- Create `test_displacement_uncertainty()` function that runs these tests and prints results
- Keep minimal verification in notebook: just call the test function once

**Action items:**
- [ ] Add test functions to DisplacementUncertainty module
- [ ] Replace TEST 2-3 cells with single call to test function
- [ ] Keep tests runnable for future validation

### Issue: TEST 1 uses VN-aligned data, doesn't belong in displacement_uncertainty
**Current:** TEST 1 estimates σ_pitch/σ_roll from VN data in the turbulence notebook

**Proposal:**
- Move TEST 1 to `lidar_vn_sync_workbench.ipynb` where VN synchronization is developed/tested
- This keeps VN-specific code with VN processing
- Turbulence notebook only needs the final σ estimates, not the data collection

**Action items:**
- [ ] Move TEST 1 cell to lidar_vn_sync_workbench.ipynb
- [ ] Document expected σ_pitch/σ_roll values for different data quality scenarios
- [ ] Turbulence notebook uses estimate_pitch_roll_uncertainty() directly in production loop

## Module Organization

### Issue: Methods duplicated between DopplerTurbulence.jl and notebook

**Current state:**
- `displacements()`, `rhopair()`, `D2_rho_stare()` etc. defined in DopplerTurbulence.jl
- Same/similar methods redefined in lidar_turbulence_production.ipynb
- Notebook versions likely newer/better

**Analysis needed:**
- [ ] Compare `displacements()` in module vs notebook
- [ ] Compare `D2_rho_stare()` in module vs notebook  
- [ ] Compare `equal_bin()`, `binavgvar()` availability
- [ ] Identify which version is canonical

**Proposal:**
1. **Consolidate to DopplerTurbulence.jl:**
   - Keep newest/working version of each method
   - Export from module
   - Import and use in notebook (don't redefine)

2. **Benefits:**
   - Single source of truth
   - Version control easier
   - Reusable across notebooks
   - Compile-time optimization

**Action items:**
- [ ] Audit method versions (module vs notebook)
- [ ] Move canonical versions to DopplerTurbulence.jl
- [ ] Update module exports
- [ ] Remove redefinitions from notebook
- [ ] Test that notebook still works with imported methods

### Issue: Check other dependents of DopplerTurbulence.jl

**Question:** What other files import DopplerTurbulence?

**Action items:**
- [ ] `grep -r "using.*DopplerTurbulence" .` to find imports
- [ ] Check if changes break other notebooks/scripts
- [ ] Update imports if module exports change

## Import Organization

### Issue: DisplacementUncertainty imported mid-notebook

**Current:**
- `include("displacement_uncertainty.jl")` in TEST 1 cell
- Not available in preamble

**Proposal:**
- Move to notebook preamble (with other includes/imports)
- Standard structure: includes at top, then using statements

**Two options:**

**Option A: Direct import (current approach)**
```julia
# Preamble
include("displacement_uncertainty.jl")
using .DisplacementUncertainty
```

**Option B: Through DopplerTurbulence (cleaner)**
```julia
# In DopplerTurbulence.jl:
include("displacement_uncertainty.jl")
using .DisplacementUncertainty
export displacement_variance_pitch_roll, propagate_rho_uncertainty, ...

# In notebook preamble:
using DopplerTurbulence  # gets DisplacementUncertainty exports too
```

**Recommendation:** Option B
- Single module to import
- DisplacementUncertainty is conceptually part of turbulence analysis
- Cleaner notebook preamble

**Action items:**
- [ ] Add DisplacementUncertainty to DopplerTurbulence.jl exports
- [ ] Move include/using to notebook preamble
- [ ] Remove from TEST cells

## Implementation Order

1. **Now:** Continue with sigma_pitch/roll implementation (don't break flow)
2. **After implementation works:** Do cleanup
3. **Test after each cleanup step:** Verify notebook still runs

## Additional Issues Found During Testing

### Overlapping method definitions (non-exported)
**Issue:** DopplerTurbulence.jl and notebook cell 3 both define some methods (e.g., utility functions, data structures) but module versions are not exported. Works at runtime but confusing for programmers.

**Action items:**
- [ ] Audit which methods overlap between DopplerTurbulence.jl and notebook cells
- [ ] Decide canonical location for each (module vs notebook)
- [ ] Remove duplicates or clearly document why both exist
- [ ] Consider: should utility functions be in a separate Utils module?

## Success Criteria

- [ ] All tests pass after cleanup
- [ ] No code duplication between module and notebook
- [ ] Clear module boundaries (DopplerTurbulence owns displacement/structure function code)
- [ ] Imports organized at top of notebook
- [ ] DisplacementUncertainty accessible through DopplerTurbulence
- [ ] No confusing overlapping method definitions
