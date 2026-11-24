# Graph Report - ASTRAL2024  (2026-11-24)

## Corpus Check
- 50 files · ~57,113 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 750 nodes · 1011 edges · 71 communities (52 shown, 18 thin omitted)
- Extraction: 98% EXTRACTED · 2% INFERRED · 0% AMBIGUOUS · INFERRED: 23 edges (avg confidence: 0.87)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `b047ef58`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- read_lidar
- LidarVNSync
- lidar_turbulence_cleanup.jl
- NoaaDas
- timing_lidar
- DopplerTurbulence
- Turbulence processing
- read_stare_chunk
- read_lidar Offset Fitting
- extract_sync_window
- vectornav2jld.jl
- time_tests.jl
- compute_hmix_epsilon.jl
- read_lidar
- read_lidar
- lidar_index
- chunk_failure_scan.jl
- Cleanup Plan After Pitch/Roll Uncertainty Implementation
- foo.jl
- coarse_and_fine_lag
- D2_rho_stare
- small_production_script.jl
- save_lidar_index.jl
- epsi_jld2nc.jl
- fft_derivative_test.jl
- Timing Offset Equations & Prior Reset Strategy
- save_lidar_dt.jl
- rain_flag_scan.jl
- Executive Summary: Timing Offset Double-Count Bug Fix
- save_stare_hourly_netcdf.jl
- Agent Init Table
- graphify
- iterative_heave_rain.jl
- lidar_uv.jl
- chunk_lidar_datetimes
- setup_sync_context_nc
- 18s GPS-UTC Leapsecond Offset (May-June 2024)
- graphify
- save_stare_mdv_sync_pass2.jl
- Alex
- FnISO and Cyclone Remal (2024)
- Leah
- Caveman Mode (AGENTS.md)
- Agent Routing & Workflow
- Persona Routing: Orchestrator to Planner to Editor to Tester
- Editor
- Orchestrator
- extend_pycall_missing.jl
- Planner
- EKAMSAT ABL Working Group Meeting (2026-04-09)
- CLAUDE.md - Project Instructions
- process_sync_data
- Tester
- Leg 1
- Production Run Additions
- f
- offset_cov
- ABLwh20260409.md
- offset_cov
- plot_tke_dissipation_may18_june09.jl
- avg_good
- filt
- load_se_epsilon_jld
- close
- check-graphify-deps.sh
- run-in-lidar-env.sh
- trigs
- recompute_mdv_fallback.jl
- ShipPosmv
- test_mdv_fallback.jl

## God Nodes (most connected - your core abstractions)
1. `LidarVNSync` - 76 edges
2. `DopplerTurbulence` - 50 edges
3. `timing_lidar` - 29 edges
4. `read_lidar` - 23 edges
5. `chunks` - 19 edges
6. `extract_sync_window()` - 16 edges
7. `refine_offset_20hz()` - 15 edges
8. `stare` - 14 edges
9. `NoaaDas` - 14 edges
10. `coarse_and_fine_lag()` - 13 edges

## Surprising Connections (you probably didn't know these)
- `Test Organization` --references--> `displacement_variance_pitch_roll()`  [EXTRACTED]
  cleanup_after_sigma_pitch_roll.md → displacement_uncertainty.jl
- `Test Organization` --references--> `propagate_rho_uncertainty()`  [EXTRACTED]
  cleanup_after_sigma_pitch_roll.md → displacement_uncertainty.jl
- `Test Organization` --references--> `estimate_pitch_roll_uncertainty()`  [EXTRACTED]
  cleanup_after_sigma_pitch_roll.md → displacement_uncertainty.jl
- `Timing Offset Composition Equations` --references--> `chunk_lidar_datetimes()`  [EXTRACTED]
  .github/plans/OFFSET_EQUATIONS_AND_RESET_STRATEGY.md → lidar_vn_sync.jl
- `Timing Offset Composition Equations` --references--> `fft_xcorr_lag()`  [EXTRACTED]
  .github/plans/OFFSET_EQUATIONS_AND_RESET_STRATEGY.md → lidar_vn_sync.jl

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Four-Persona Agent Workflow (Orchestrator/Planner/Editor/Tester)** — github_agents_orchestrator_agent_orchestrator, github_agents_planner_agent_planner, github_agents_editor_agent_editor, github_agents_tester_agent_tester, agents_persona_routing [EXTRACTED 1.00]
- **Timing Offset Double-Count Bug: Diagnosis, Equations, Fix, and Failure-Mode Rationale** — github_plans_executive_summary_offset_fix_double_count_bug, github_plans_executive_summary_offset_fix_offset_composition_fix, github_plans_offset_equations_and_reset_strategy_offset_equations, github_plans_offset_equations_and_reset_strategy_prior_reset_strategy, todo_sync_failure_modes [INFERRED 0.85]
- **Structure-Function Dissipation Retrieval Pipeline** — readme_structure_function_fit, readme_epsilon_retrieval, production_add_qc_fit_stats_onepass, production_add_qc_se_epsilon_propagation, production_add_qc_epsilon_sentinels [INFERRED 0.85]

## Communities (71 total, 18 thin omitted)

### Community 0 - "read_lidar"
Cohesion: 0.08
Nodes (32): all_chunks(), all_gaps(), all_start_end_indices(), chunken(), chunkst(), compute_mdv_snr_mean(), get_all_file_start_end_idxs(), get_daily_meanuv() (+24 more)

### Community 1 - "LidarVNSync"
Cohesion: 0.07
Nodes (22): fcn3(), chunks, Dates, DSP, FFTW, Interpolations, JLD2, NCDatasets (+14 more)

### Community 2 - "lidar_turbulence_cleanup.jl"
Cohesion: 0.06
Nodes (18): chunks, Dates, DSP, FFTW, Interpolations, JLD2, NCDatasets, NoaaDas (+10 more)

### Community 3 - "NoaaDas"
Cohesion: 0.16
Nodes (27): das_dict(), DasFps, DasGps, DasScs, declat(), declon(), get_das_filenames(), get_das_pathfiles() (+19 more)

### Community 4 - "timing_lidar"
Cohesion: 0.10
Nodes (21): ba(), cross_correlation(), dtregress(), DVM(), find_lags_iterative(), gps2dt(), gps2utc(), gpstime2gpsvndt() (+13 more)

### Community 5 - "DopplerTurbulence"
Cohesion: 0.06
Nodes (39): Module Consolidation into DopplerTurbulence.jl, DisplacementUncertainty, allcross(), anom(), D2_rho_small_disp(), D2_rho_stare(), D2_rho_stare_qc(), displacements() (+31 more)

### Community 6 - "Turbulence processing"
Cohesion: 0.09
Nodes (25): Devmi, Entrainment and Cold Pools from Scalar Budgets, Jayesh, Lidar Backscatter BL Height (Haar Wavelet Gradient Detection), Mixing Depth from Lidar (Next Step), Simon, Epsilon Conversion epsilon(A), C2ll=2.0 (QC doc), Epsilon Sentinel Coding (QC doc) (+17 more)

### Community 7 - "read_stare_chunk"
Cohesion: 0.29
Nodes (8): binavg(), equal_bin(), f(), get_mdv(), indavg(), pcolor_lidar_stare(), read_stare_chunk(), remove_mdv()

### Community 8 - "read_lidar Offset Fitting"
Cohesion: 0.16
Nodes (12): chunks, dt_to_chunkind(), findindices(), finelagcov(), fit_offset(), get_mean_uv_chunk!(), goodcov(), indavg() (+4 more)

### Community 9 - "extract_sync_window"
Cohesion: 0.23
Nodes (15): diag_array(), diagnostic_example_chunks(), diagnostic_single_window(), dt_seconds(), ensure_chunk_loaded!(), extract_sync_window(), finite_overlap_corr(), init_stream_state() (+7 more)

### Community 10 - "vectornav2jld.jl"
Cohesion: 0.13
Nodes (9): ShipPosmv, Dates, Interpolations, JLD2, NoaaDas, Pkg, read_vecnav, Revise (+1 more)

### Community 11 - "time_tests.jl"
Cohesion: 0.17
Nodes (8): chunken(), chunkst(), Dates, read_lidar, Revise, nextchunki(), prevchunki(), thisj()

### Community 12 - "compute_hmix_epsilon.jl"
Cohesion: 0.24
Nodes (12): compute_hmix_for_matrix(), extract_epsilon_and_meta(), hmix_first_crossing(), is_valid_eps(), Dates, JLD2, Pkg, Statistics (+4 more)

### Community 13 - "read_lidar"
Cohesion: 0.21
Nodes (8): get_daily_meanuv(), Dates, Interpolations, JLD2, NCDatasets, read_lidar, read_streamlinexr_head(), read_streamlinexr_stare!()

### Community 14 - "read_lidar"
Cohesion: 0.21
Nodes (8): get_daily_meanuv(), Dates, Interpolations, JLD2, NCDatasets, read_lidar, read_streamlinexr_head(), read_streamlinexr_stare!()

### Community 15 - "lidar_index"
Cohesion: 0.20
Nodes (5): build_lidar_index(), Dates, JLD2, lidar_index, LidarIndex

### Community 16 - "chunk_failure_scan.jl"
Cohesion: 0.20
Nodes (7): Dates, FFTW, LidarVNSync, NCDatasets, Pkg, Printf, Statistics

### Community 17 - "Cleanup Plan After Pitch/Roll Uncertainty Implementation"
Cohesion: 0.10
Nodes (18): Additional Issues Found During Testing, Cleanup Plan After Pitch/Roll Uncertainty Implementation, Implementation Order, Import Organization, Issue: Check other dependents of DopplerTurbulence.jl, Issue: DisplacementUncertainty imported mid-notebook, Issue: Methods duplicated between DopplerTurbulence.jl and notebook, Issue: TEST 1 uses VN-aligned data, doesn't belong in displacement_uncertainty (+10 more)

### Community 18 - "foo.jl"
Cohesion: 0.31
Nodes (7): estimate_heave_ts(), iterative_despike_rain(), Statistics, nanmedian(), rain_mask_vertical(), running_median(), subtract_rain_layers!()

### Community 19 - "coarse_and_fine_lag"
Cohesion: 0.25
Nodes (11): analytic_envelope_fft(), coarse_and_fine_lag(), cosine_edge_mask(), detrend_center(), fft_bandpass(), fft_xcorr_lag(), fill_short_nan_gaps(), finite_overlap() (+3 more)

### Community 20 - "D2_rho_stare"
Cohesion: 0.29
Nodes (7): allcross(), anom(), D2_rho_stare(), displacements(), lidarindices(), rng(), uniquepairs()

### Community 21 - "small_production_script.jl"
Cohesion: 0.05
Nodes (25): Profile, PyCall, chunks, Dates, Distributions, DopplerTurbulence, DSP, FFTW (+17 more)

### Community 22 - "save_lidar_index.jl"
Cohesion: 0.25
Nodes (6): JLD2, lidar_index, NCDatasets, Pkg, read_lidar, Revise

### Community 23 - "epsi_jld2nc.jl"
Cohesion: 0.29
Nodes (4): Dates, JLD2, NCDatasets, Pkg

### Community 24 - "fft_derivative_test.jl"
Cohesion: 0.33
Nodes (6): f_taper(), FFTW, Pkg, Revise, taper(), Plots

### Community 25 - "Timing Offset Equations & Prior Reset Strategy"
Cohesion: 0.10
Nodes (20): Change 1: `run_sequential_offsets()` signature, Change 2: `prior_from_history()` logic, Change 3: Track file index in `run_sequential_offsets()`, Change 4: Fix the double-count bug in `process_sync_data()`, Correct Equations (What We Want), Current Implementation (Broken), Equation 1: Lidar time synchronization, Equation 2: 1 Hz coarse + fine lag stage (+12 more)

### Community 26 - "save_lidar_dt.jl"
Cohesion: 0.29
Nodes (6): Dates, JLD2, NCDatasets, Pkg, read_lidar, Revise

### Community 27 - "rain_flag_scan.jl"
Cohesion: 0.33
Nodes (4): FFTW, LinearAlgebra, Statistics, ImageFiltering

### Community 28 - "Executive Summary: Timing Offset Double-Count Bug Fix"
Cohesion: 0.13
Nodes (14): After (New Code):, Before (Old Code):, Code Changes, Configurable Prior Reset (if cascade failures are a problem), Executive Summary: Timing Offset Double-Count Bug Fix, Files Updated, How to Verify, Next Steps (Optional, Not Critical) (+6 more)

### Community 29 - "save_stare_hourly_netcdf.jl"
Cohesion: 0.33
Nodes (4): Dates, Pkg, Printf, read_lidar

### Community 30 - "Agent Init Table"
Cohesion: 1.00
Nodes (5): Agent Init Table, Editor Agent, Orchestrator Agent, Planner Agent, Tester Agent

### Community 31 - "graphify"
Cohesion: 0.40
Nodes (5): graphify Mandatory Resource Rules (AGENTS.md), AGENTS.md - Agent Routing & Workflow, graphify Rules (CLAUDE.md), Other Agent Contexts Note (legacy docs deprecated for Claude Code), graphify

### Community 32 - "iterative_heave_rain.jl"
Cohesion: 0.60
Nodes (4): infer_heave_from_stare(), Statistics, local_vertical_stats(), weighted_median()

### Community 33 - "lidar_uv.jl"
Cohesion: 0.40
Nodes (4): NCDatasets, Pkg, PythonPlot, Statistics

### Community 34 - "chunk_lidar_datetimes"
Cohesion: 0.67
Nodes (4): check_chunk_alignment_contract(), chunk_lidar_datetimes(), ensure_chunk_loaded_nc!(), load_chunk_for_dissipation()

### Community 35 - "setup_sync_context_nc"
Cohesion: 0.50
Nodes (5): init_periodic_beams(), load_lidar_indices_and_files(), setup_sync_context(), setup_sync_context_nc(), write_mdv_sync_pass2!()

### Community 36 - "18s GPS-UTC Leapsecond Offset (May-June 2024)"
Cohesion: 0.50
Nodes (4): 18s GPS-UTC Leapsecond Offset (May-June 2024), Leg 1 VN-POSMV Clock Sync (13-hr precession), Leg 2 VN-POSMV Sync (0.6s lead, discretization artifact), Motion Timing Report (VectorNav/POSMV)

### Community 43 - "Agent Routing & Workflow"
Cohesion: 0.07
Nodes (27): Agent Descriptions & When LLMs Must Use Them, Agent Init, Agent Routing & Workflow, Caveman Mode (Active), Editor, Environment, graphify (Mandatory Resource), Key data structures (+19 more)

### Community 45 - "Editor"
Cohesion: 0.33
Nodes (5): Code Style Context, Core Responsibilities, Editor, Output Format (Caveman Clause), When to Use This Agent

### Community 46 - "Orchestrator"
Cohesion: 0.33
Nodes (5): Core Responsibilities, Implementation Style, Orchestrator, Output Format (Caveman Clause), When to Use This Agent

### Community 48 - "Planner"
Cohesion: 0.33
Nodes (5): Code Style Context, Core Responsibilities, Output Format (Caveman Clause), Planner, When to Use This Agent

### Community 51 - "process_sync_data"
Cohesion: 0.17
Nodes (16): LidarVNSync Core API Conventions (AGENTS.md), LidarVNSync Core API Conventions (CLAUDE.md), Prior-Offset Double-Count Bug, Offset Composition Fix (single-residual rule), Timing Offset Composition Equations, Prior Reset Strategy (carry_forward / reset_each_file / reset_on_bad_quality), append_nan_offset_log(), backward_jump_robustness() (+8 more)

### Community 52 - "Tester"
Cohesion: 0.33
Nodes (5): Best Practices, Core Responsibilities, Output Format (Caveman Clause), Tester, When to Use This Agent

### Community 53 - "Leg 1"
Cohesion: 0.33
Nodes (5): A Report gleaned from experiments in [timing_lidar.jl](timing_lidar.jl) and [lidar_turbulence.jl](lidar_turbulence.jl), Caption: (+) Positive lag means Vn start index is shifted forward, i.e., Vn signal lags POSMV. (Data in this shifted index is (moved backward and) aligned with the unshifted POSMV start.) (-) Negative lag means negative index of Vn lines up with start of POSMV epoch, i.e., Vn signal leads POSMV., Leg 1, Leg 2, part 1, Timing of motion data for the Halo Photonics lidar from VectorNav and POSMV

### Community 54 - "Production Run Additions"
Cohesion: 0.33
Nodes (5): Efficiency changes, Output variables saved, Postproduction utility, Principles of calculation, Production Run Additions

### Community 55 - "f"
Cohesion: 0.40
Nodes (5): binavgvar(), f(), get_mdv(), pcolor_lidar_stare(), remove_mdv()

### Community 56 - "offset_cov"
Cohesion: 0.40
Nodes (5): good(), offset_cov(), offset_range_covs(), offset_subset(), sync_offset()

### Community 57 - "ABLwh20260409.md"
Cohesion: 0.50
Nodes (3): Agenda, EKAMSAT ABL Working Group, Lidar turbulence analysis - Simon

### Community 58 - "offset_cov"
Cohesion: 0.33
Nodes (6): findindices(), good(), offset_cov(), offset_range_covs(), offset_subset(), sync_offset()

### Community 59 - "plot_tke_dissipation_may18_june09.jl"
Cohesion: 0.33
Nodes (5): Dates, JLD2, Pkg, PythonPlot, Statistics

### Community 68 - "recompute_mdv_fallback.jl"
Cohesion: 0.21
Nodes (17): aligned_vn(), blend(), chunk_winds(), daily_wind!(), fit_epsilon(), good(), interpolate_at(), Dates (+9 more)

### Community 69 - "ShipPosmv"
Cohesion: 0.39
Nodes (9): Base.Iterators, get_nav_file(), get_posmv_file(), itr_expand(), read_gyro_data(), read_gyro_dict(), read_pashr_data(), read_pashr_dict() (+1 more)

## Knowledge Gaps
- **264 isolated node(s):** `check-graphify-deps.sh script`, `run-in-lidar-env.sh script`, `bash`, `Revise`, `Pkg` (+259 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 401 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **18 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `LidarVNSync` connect `LidarVNSync` to `chunk_lidar_datetimes`, `setup_sync_context_nc`, `extract_sync_window`, `process_sync_data`, `coarse_and_fine_lag`?**
  _High betweenness centrality (0.016) - this node is a cross-community bridge._
- **Why does `Module Consolidation into DopplerTurbulence.jl` connect `DopplerTurbulence` to `Cleanup Plan After Pitch/Roll Uncertainty Implementation`, `Turbulence processing`?**
  _High betweenness centrality (0.011) - this node is a cross-community bridge._
- **What connects `check-graphify-deps.sh script`, `run-in-lidar-env.sh script`, `bash` to the rest of the system?**
  _264 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `read_lidar` be split into smaller, more focused modules?**
  _Cohesion score 0.08170731707317073 - nodes in this community are weakly interconnected._
- **Should `LidarVNSync` be split into smaller, more focused modules?**
  _Cohesion score 0.06818181818181818 - nodes in this community are weakly interconnected._
- **Should `lidar_turbulence_cleanup.jl` be split into smaller, more focused modules?**
  _Cohesion score 0.058823529411764705 - nodes in this community are weakly interconnected._
- **Should `timing_lidar` be split into smaller, more focused modules?**
  _Cohesion score 0.10344827586206896 - nodes in this community are weakly interconnected._