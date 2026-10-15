# Graph Report - ASTRAL2024  (2026-10-15)

## Corpus Check
- 47 files · ~54,613 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 700 nodes · 947 edges · 65 communities (49 shown, 16 thin omitted)
- Extraction: 98% EXTRACTED · 2% INFERRED · 0% AMBIGUOUS · INFERRED: 23 edges (avg confidence: 0.87)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `de3a05d5`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- [[_COMMUNITY_read_lidar|read_lidar]]
- [[_COMMUNITY_LidarVNSync|LidarVNSync]]
- [[_COMMUNITY_lidar_turbulence_cleanup.jl|lidar_turbulence_cleanup.jl]]
- [[_COMMUNITY_NoaaDas|NoaaDas]]
- [[_COMMUNITY_timing_lidar|timing_lidar]]
- [[_COMMUNITY_DopplerTurbulence|DopplerTurbulence]]
- [[_COMMUNITY_Saved JLD Output Variables (QC doc)|Saved JLD Output Variables (QC doc)]]
- [[_COMMUNITY_run_sequential_offsets|run_sequential_offsets]]
- [[_COMMUNITY_read_lidar Offset Fitting|read_lidar Offset Fitting]]
- [[_COMMUNITY_extract_sync_window|extract_sync_window]]
- [[_COMMUNITY_vectornav2jld.jl|vectornav2jld.jl]]
- [[_COMMUNITY_time_tests.jl|time_tests.jl]]
- [[_COMMUNITY_compute_hmix_epsilon.jl|compute_hmix_epsilon.jl]]
- [[_COMMUNITY_read_lidar|read_lidar]]
- [[_COMMUNITY_read_lidar|read_lidar]]
- [[_COMMUNITY_lidar_index|lidar_index]]
- [[_COMMUNITY_chunk_failure_scan.jl|chunk_failure_scan.jl]]
- [[_COMMUNITY_DisplacementUncertainty|DisplacementUncertainty]]
- [[_COMMUNITY_foo.jl|foo.jl]]
- [[_COMMUNITY_coarse_and_fine_lag|coarse_and_fine_lag]]
- [[_COMMUNITY_ShipPosmv|ShipPosmv]]
- [[_COMMUNITY_read_stare_chunk|read_stare_chunk]]
- [[_COMMUNITY_save_lidar_index.jl|save_lidar_index.jl]]
- [[_COMMUNITY_epsi_jld2nc.jl|epsi_jld2nc.jl]]
- [[_COMMUNITY_fft_derivative_test.jl|fft_derivative_test.jl]]
- [[_COMMUNITY_D2_rho_stare|D2_rho_stare]]
- [[_COMMUNITY_save_lidar_dt.jl|save_lidar_dt.jl]]
- [[_COMMUNITY_rain_flag_scan.jl|rain_flag_scan.jl]]
- [[_COMMUNITY_offset_cov|offset_cov]]
- [[_COMMUNITY_save_stare_hourly_netcdf.jl|save_stare_hourly_netcdf.jl]]
- [[_COMMUNITY_Agent Init Table|Agent Init Table]]
- [[_COMMUNITY_Graphify Workflow (Copilot Instructions)|Graphify Workflow (Copilot Instructions)]]
- [[_COMMUNITY_iterative_heave_rain.jl|iterative_heave_rain.jl]]
- [[_COMMUNITY_lidar_uv.jl|lidar_uv.jl]]
- [[_COMMUNITY_chunk_lidar_datetimes|chunk_lidar_datetimes]]
- [[_COMMUNITY_init_stream_state|init_stream_state]]
- [[_COMMUNITY_18s GPS-UTC Leapsecond Offset (May-June 2024)|18s GPS-UTC Leapsecond Offset (May-June 2024)]]
- [[_COMMUNITY_graphify|graphify]]
- [[_COMMUNITY_save_stare_mdv_sync_pass2.jl|save_stare_mdv_sync_pass2.jl]]
- [[_COMMUNITY_Alex|Alex]]
- [[_COMMUNITY_FnISO and Cyclone Remal (2024)|FnISO and Cyclone Remal (2024)]]
- [[_COMMUNITY_Leah|Leah]]
- [[_COMMUNITY_Caveman Mode (AGENTS.md)|Caveman Mode (AGENTS.md)]]
- [[_COMMUNITY_lidar micromambaconda Environment (AGENTS.md)|lidar micromamba/conda Environment (AGENTS.md)]]
- [[_COMMUNITY_Persona Routing Orchestrator to Planner to Editor to Tester|Persona Routing: Orchestrator to Planner to Editor to Tester]]
- [[_COMMUNITY_Graphify Deps Check Hook|Graphify Deps Check Hook]]
- [[_COMMUNITY_Lidar Env Wrapper Hook|Lidar Env Wrapper Hook]]
- [[_COMMUNITY_extend_pycall_missing.jl|extend_pycall_missing.jl]]
- [[_COMMUNITY_trigs|trigs]]
- [[_COMMUNITY_EKAMSAT ABL Working Group Meeting (2026-04-09)|EKAMSAT ABL Working Group Meeting (2026-04-09)]]
- [[_COMMUNITY_CLAUDE.md - Project Instructions|CLAUDE.md - Project Instructions]]
- [[_COMMUNITY_process_sync_data|process_sync_data]]
- [[_COMMUNITY_Community 52|Community 52]]
- [[_COMMUNITY_Community 53|Community 53]]
- [[_COMMUNITY_Community 54|Community 54]]
- [[_COMMUNITY_Community 55|Community 55]]
- [[_COMMUNITY_Community 56|Community 56]]
- [[_COMMUNITY_Community 57|Community 57]]
- [[_COMMUNITY_Community 58|Community 58]]
- [[_COMMUNITY_Community 59|Community 59]]
- [[_COMMUNITY_Community 60|Community 60]]
- [[_COMMUNITY_Community 61|Community 61]]
- [[_COMMUNITY_Community 62|Community 62]]
- [[_COMMUNITY_Community 63|Community 63]]

## God Nodes (most connected - your core abstractions)
1. `LidarVNSync` - 74 edges
2. `DopplerTurbulence` - 45 edges
3. `timing_lidar` - 29 edges
4. `read_lidar` - 23 edges
5. `chunks` - 19 edges
6. `extract_sync_window()` - 16 edges
7. `refine_offset_20hz()` - 15 edges
8. `stare` - 14 edges
9. `NoaaDas` - 14 edges
10. `coarse_and_fine_lag()` - 13 edges

## Surprising Connections (you probably didn't know these)
- `Sync Failure Modes & Prior Handling` --semantically_similar_to--> `Prior-Offset Double-Count Bug`  [INFERRED] [semantically similar]
  TODO.md → .github/plans/EXECUTIVE_SUMMARY_OFFSET_FIX.md
- `Test Organization` --references--> `displacement_variance_pitch_roll()`  [EXTRACTED]
  cleanup_after_sigma_pitch_roll.md → displacement_uncertainty.jl
- `Test Organization` --references--> `propagate_rho_uncertainty()`  [EXTRACTED]
  cleanup_after_sigma_pitch_roll.md → displacement_uncertainty.jl
- `Test Organization` --references--> `estimate_pitch_roll_uncertainty()`  [EXTRACTED]
  cleanup_after_sigma_pitch_roll.md → displacement_uncertainty.jl
- `Timing Offset Composition Equations` --references--> `chunk_lidar_datetimes()`  [EXTRACTED]
  .github/plans/OFFSET_EQUATIONS_AND_RESET_STRATEGY.md → lidar_vn_sync.jl

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Four-Persona Agent Workflow (Orchestrator/Planner/Editor/Tester)** — github_agents_orchestrator_agent_orchestrator, github_agents_planner_agent_planner, github_agents_editor_agent_editor, github_agents_tester_agent_tester, agents_persona_routing [EXTRACTED 1.00]
- **Timing Offset Double-Count Bug: Diagnosis, Equations, Fix, and Failure-Mode Rationale** — github_plans_executive_summary_offset_fix_double_count_bug, github_plans_executive_summary_offset_fix_offset_composition_fix, github_plans_offset_equations_and_reset_strategy_offset_equations, github_plans_offset_equations_and_reset_strategy_prior_reset_strategy, todo_sync_failure_modes [INFERRED 0.85]
- **Structure-Function Dissipation Retrieval Pipeline** — readme_structure_function_fit, readme_epsilon_retrieval, production_add_qc_fit_stats_onepass, production_add_qc_se_epsilon_propagation, production_add_qc_epsilon_sentinels [INFERRED 0.85]

## Communities (65 total, 16 thin omitted)

### Community 0 - "read_lidar"
Cohesion: 0.08
Nodes (32): all_chunks(), all_gaps(), all_start_end_indices(), chunken(), chunkst(), compute_mdv_snr_mean(), get_all_file_start_end_idxs(), get_daily_meanuv() (+24 more)

### Community 1 - "LidarVNSync"
Cohesion: 0.07
Nodes (20): fcn3(), Dates, DSP, FFTW, Interpolations, JLD2, NCDatasets, Printf (+12 more)

### Community 2 - "lidar_turbulence_cleanup.jl"
Cohesion: 0.05
Nodes (38): allcross(), anom(), binavg(), D2_rho_stare(), displacements(), equal_bin(), f(), findindices() (+30 more)

### Community 3 - "NoaaDas"
Cohesion: 0.15
Nodes (27): das_dict, das_dict(), DasFps, DasGps, DasScs, declat(), declon(), get_das_filenames() (+19 more)

### Community 4 - "timing_lidar"
Cohesion: 0.10
Nodes (21): ba(), cross_correlation(), dtregress(), DVM(), find_lags_iterative(), gps2dt(), gps2utc(), gpstime2gpsvndt() (+13 more)

### Community 5 - "DopplerTurbulence"
Cohesion: 0.07
Nodes (34): Module Consolidation into DopplerTurbulence.jl, allcross(), anom(), D2_rho_small_disp(), D2_rho_stare(), D2_rho_stare_qc(), displacements(), DopplerTurbulence (+26 more)

### Community 6 - "Saved JLD Output Variables (QC doc)"
Cohesion: 0.09
Nodes (25): Devmi, Entrainment and Cold Pools from Scalar Budgets, Jayesh, Lidar Backscatter BL Height (Haar Wavelet Gradient Detection), Mixing Depth from Lidar (Next Step), Simon, Epsilon Conversion epsilon(A), C2ll=2.0 (QC doc), Epsilon Sentinel Coding (QC doc) (+17 more)

### Community 7 - "run_sequential_offsets"
Cohesion: 0.29
Nodes (10): Sync-quality gates, Sync-quality gates, Prior Reset Strategy (carry_forward / reset_each_file / reset_on_bad_quality), backward_jump_robustness(), is_jump_candidate(), prior_from_history(), robust_center(), run_sequential_offsets() (+2 more)

### Community 8 - "read_lidar Offset Fitting"
Cohesion: 0.16
Nodes (12): chunks, dt_to_chunkind(), findindices(), finelagcov(), fit_offset(), get_mean_uv_chunk!(), goodcov(), indavg() (+4 more)

### Community 9 - "extract_sync_window"
Cohesion: 0.18
Nodes (18): diag_array(), diagnostic_example_chunks(), diagnostic_single_window(), dt_seconds(), ensure_chunk_loaded!(), extract_sync_window(), fill_short_nan_gaps(), finite_overlap() (+10 more)

### Community 10 - "vectornav2jld.jl"
Cohesion: 0.17
Nodes (6): Dates, Interpolations, JLD2, Pkg, read_vecnav_dict, Revise

### Community 11 - "time_tests.jl"
Cohesion: 0.19
Nodes (7): chunken(), chunkst(), Dates, Revise, nextchunki(), prevchunki(), thisj()

### Community 12 - "compute_hmix_epsilon.jl"
Cohesion: 0.24
Nodes (12): compute_hmix_for_matrix(), extract_epsilon_and_meta(), hmix_first_crossing(), is_valid_eps(), Dates, JLD2, Pkg, Statistics (+4 more)

### Community 13 - "read_lidar"
Cohesion: 0.21
Nodes (8): get_daily_meanuv(), read_lidar, read_streamlinexr_head(), read_streamlinexr_stare!(), Dates, Interpolations, JLD2, NCDatasets

### Community 14 - "read_lidar"
Cohesion: 0.21
Nodes (8): get_daily_meanuv(), read_lidar, read_streamlinexr_head(), read_streamlinexr_stare!(), Dates, Interpolations, JLD2, NCDatasets

### Community 15 - "lidar_index"
Cohesion: 0.20
Nodes (5): build_lidar_index(), Dates, JLD2, lidar_index, LidarIndex

### Community 16 - "chunk_failure_scan.jl"
Cohesion: 0.22
Nodes (6): Dates, FFTW, NCDatasets, Pkg, Printf, Statistics

### Community 17 - "DisplacementUncertainty"
Cohesion: 0.10
Nodes (18): Additional Issues Found During Testing, Cleanup Plan After Pitch/Roll Uncertainty Implementation, Implementation Order, Import Organization, Issue: Check other dependents of DopplerTurbulence.jl, Issue: DisplacementUncertainty imported mid-notebook, Issue: Methods duplicated between DopplerTurbulence.jl and notebook, Issue: TEST 1 uses VN-aligned data, doesn't belong in displacement_uncertainty (+10 more)

### Community 18 - "foo.jl"
Cohesion: 0.31
Nodes (7): estimate_heave_ts(), iterative_despike_rain(), Statistics, nanmedian(), rain_mask_vertical(), running_median(), subtract_rain_layers!()

### Community 19 - "coarse_and_fine_lag"
Cohesion: 0.31
Nodes (10): Timing Offset Composition Equations, analytic_envelope_fft(), coarse_and_fine_lag(), cosine_edge_mask(), detrend_center(), fft_bandpass(), fft_xcorr_lag(), iterative_coarse_lag() (+2 more)

### Community 20 - "ShipPosmv"
Cohesion: 0.39
Nodes (9): flatten, get_nav_file(), get_posmv_file(), itr_expand(), read_gyro_data(), read_gyro_dict(), read_pashr_data(), read_pashr_dict() (+1 more)

### Community 21 - "read_stare_chunk"
Cohesion: 0.05
Nodes (21): Distributions, Profile, PyCall, cat_dicts, Dates, DSP, FFTW, Interpolations (+13 more)

### Community 22 - "save_lidar_index.jl"
Cohesion: 0.33
Nodes (4): JLD2, NCDatasets, Pkg, Revise

### Community 23 - "epsi_jld2nc.jl"
Cohesion: 0.29
Nodes (4): Dates, JLD2, NCDatasets, Pkg

### Community 24 - "fft_derivative_test.jl"
Cohesion: 0.33
Nodes (6): f_taper(), FFTW, Pkg, Revise, taper(), Plots

### Community 25 - "D2_rho_stare"
Cohesion: 0.10
Nodes (20): Change 1: `run_sequential_offsets()` signature, Change 2: `prior_from_history()` logic, Change 3: Track file index in `run_sequential_offsets()`, Change 4: Fix the double-count bug in `process_sync_data()`, Correct Equations (What We Want), Current Implementation (Broken), Equation 1: Lidar time synchronization, Equation 2: 1 Hz coarse + fine lag stage (+12 more)

### Community 26 - "save_lidar_dt.jl"
Cohesion: 0.33
Nodes (5): Dates, JLD2, NCDatasets, Pkg, Revise

### Community 27 - "rain_flag_scan.jl"
Cohesion: 0.33
Nodes (4): FFTW, LinearAlgebra, Statistics, ImageFiltering

### Community 28 - "offset_cov"
Cohesion: 0.13
Nodes (14): After (New Code):, Before (Old Code):, Code Changes, Configurable Prior Reset (if cascade failures are a problem), Executive Summary: Timing Offset Double-Count Bug Fix, Files Updated, How to Verify, Next Steps (Optional, Not Critical) (+6 more)

### Community 29 - "save_stare_hourly_netcdf.jl"
Cohesion: 0.33
Nodes (4): read_streamlinexr_stare, Dates, Pkg, Printf

### Community 30 - "Agent Init Table"
Cohesion: 1.00
Nodes (5): Agent Init Table, Editor Agent, Orchestrator Agent, Planner Agent, Tester Agent

### Community 31 - "Graphify Workflow (Copilot Instructions)"
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

### Community 35 - "init_stream_state"
Cohesion: 0.67
Nodes (4): init_periodic_beams(), load_lidar_indices_and_files(), setup_sync_context(), setup_sync_context_nc()

### Community 36 - "18s GPS-UTC Leapsecond Offset (May-June 2024)"
Cohesion: 0.50
Nodes (4): 18s GPS-UTC Leapsecond Offset (May-June 2024), Leg 1 VN-POSMV Clock Sync (13-hr precession), Leg 2 VN-POSMV Sync (0.6s lead, discretization artifact), Motion Timing Report (VectorNav/POSMV)

### Community 43 - "lidar micromamba/conda Environment (AGENTS.md)"
Cohesion: 0.08
Nodes (24): Agent Descriptions & When LLMs Must Use Them, Agent Init, Agent Routing & Workflow, Caveman Mode (Active), Editor, Environment, graphify (Mandatory Resource), Key data structures (+16 more)

### Community 45 - "Graphify Deps Check Hook"
Cohesion: 0.33
Nodes (5): Code Style Context, Core Responsibilities, Editor, Output Format (Caveman Clause), When to Use This Agent

### Community 46 - "Lidar Env Wrapper Hook"
Cohesion: 0.33
Nodes (5): Core Responsibilities, Implementation Style, Orchestrator, Output Format (Caveman Clause), When to Use This Agent

### Community 48 - "trigs"
Cohesion: 0.33
Nodes (5): Code Style Context, Core Responsibilities, Output Format (Caveman Clause), Planner, When to Use This Agent

### Community 51 - "process_sync_data"
Cohesion: 0.25
Nodes (8): LidarVNSync Core API Conventions (AGENTS.md), LidarVNSync Core API Conventions (CLAUDE.md), Prior-Offset Double-Count Bug, Offset Composition Fix (single-residual rule), append_nan_offset_log(), offset_or_sentinel(), offset_to_ms(), process_sync_data()

### Community 52 - "Community 52"
Cohesion: 0.33
Nodes (5): Best Practices, Core Responsibilities, Output Format (Caveman Clause), Tester, When to Use This Agent

### Community 53 - "Community 53"
Cohesion: 0.33
Nodes (5): A Report gleaned from experiments in [timing_lidar.jl](timing_lidar.jl) and [lidar_turbulence.jl](lidar_turbulence.jl), Caption: (+) Positive lag means Vn start index is shifted forward, i.e., Vn signal lags POSMV. (Data in this shifted index is (moved backward and) aligned with the unshifted POSMV start.) (-) Negative lag means negative index of Vn lines up with start of POSMV epoch, i.e., Vn signal leads POSMV., Leg 1, Leg 2, part 1, Timing of motion data for the Halo Photonics lidar from VectorNav and POSMV

### Community 54 - "Community 54"
Cohesion: 0.33
Nodes (5): Efficiency changes, Output variables saved, Postproduction utility, Principles of calculation, Production Run Additions

### Community 55 - "Community 55"
Cohesion: 0.40
Nodes (5): binavgvar(), f(), get_mdv(), pcolor_lidar_stare(), remove_mdv()

### Community 56 - "Community 56"
Cohesion: 0.40
Nodes (5): good(), offset_cov(), offset_range_covs(), offset_subset(), sync_offset()

### Community 57 - "Community 57"
Cohesion: 0.50
Nodes (3): Agenda, EKAMSAT ABL Working Group, Lidar turbulence analysis - Simon

## Knowledge Gaps
- **229 isolated node(s):** `check-graphify-deps.sh script`, `run-in-lidar-env.sh script`, `bash`, `Revise`, `Pkg` (+224 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **16 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `LidarVNSync` connect `LidarVNSync` to `chunk_lidar_datetimes`, `init_stream_state`, `run_sequential_offsets`, `extract_sync_window`, `process_sync_data`, `coarse_and_fine_lag`?**
  _High betweenness centrality (0.017) - this node is a cross-community bridge._
- **Why does `Module Consolidation into DopplerTurbulence.jl` connect `DopplerTurbulence` to `DisplacementUncertainty`, `Saved JLD Output Variables (QC doc)`?**
  _High betweenness centrality (0.011) - this node is a cross-community bridge._
- **What connects `check-graphify-deps.sh script`, `run-in-lidar-env.sh script`, `bash` to the rest of the system?**
  _231 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `read_lidar` be split into smaller, more focused modules?**
  _Cohesion score 0.08170731707317073 - nodes in this community are weakly interconnected._
- **Should `LidarVNSync` be split into smaller, more focused modules?**
  _Cohesion score 0.07311827956989247 - nodes in this community are weakly interconnected._
- **Should `lidar_turbulence_cleanup.jl` be split into smaller, more focused modules?**
  _Cohesion score 0.05241090146750524 - nodes in this community are weakly interconnected._
- **Should `timing_lidar` be split into smaller, more focused modules?**
  _Cohesion score 0.10344827586206896 - nodes in this community are weakly interconnected._