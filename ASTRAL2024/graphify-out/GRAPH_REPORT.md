# Graph Report - ASTRAL2024  (2026-09-07)

## Corpus Check
- Corpus is ~45,174 words - fits in a single context window. You may not need a graph.

## Summary
- 521 nodes · 740 edges · 51 communities (37 shown, 14 thin omitted)
- Extraction: 97% EXTRACTED · 3% INFERRED · 0% AMBIGUOUS · INFERRED: 23 edges (avg confidence: 0.87)
- Token cost: 0 input · 117,491 output

## Community Hubs (Navigation)
- read_lidar Core Chunking
- LidarVNSync Core
- Legacy Turbulence Cleanup (old/)
- NoaaDas Readers
- Timing Offset Analysis
- DopplerTurbulence Module
- ABL Working Group & Production QC
- Sync API Docs & Offset Sentinel
- read_lidar Offset Fitting
- Sync Window Extraction
- VectorNav to JLD2 Conversion
- Chunk Timing Tests
- Hmix/Epsilon Computation
- Deprecated Reader (2024-08-07)
- Deprecated Reader (buggy offsets)
- Lidar Index Building
- Chunk Failure Scan
- Pitch/Roll Uncertainty & Cleanup Plan
- Rain Despike Ideas (codeIdeas/foo.jl)
- FFT Lag/Xcorr & Offset Equations
- POSMV/Gyro Readers
- Legacy Turbulence Binning/Plot
- Save Lidar Index Script
- Epsilon JLD2 to NetCDF
- FFT Derivative Test (old/)
- Legacy D2(rho) Structure Function
- Save Lidar Datetime Script
- Rain Flag Scan Idea
- Legacy Offset Covariance Calc
- Save Hourly Stare NetCDF
- Copilot 4-Persona Agents (legacy)
- graphify Rules (Docs Cross-Ref)
- Iterative Heave/Rain Idea
- Lidar UV Script
- Chunk Alignment & Loading
- Sync Context Setup
- Motion Timing Doc (Leg1/Leg2)
- graphify MCP Config
- MDV Sync Pass2 Script
- ABL Meeting: Radiosonde Array
- ABL Meeting: Cyclone Remal
- ABL Meeting: MSE Budgets
- Caveman Mode (Docs)
- Env Activation Docs
- Persona Routing vs Native Loop
- Graphify Deps Check Hook
- Lidar Env Wrapper Hook
- PyCall Extension Script
- Legacy Trig Helpers
- ABL Meeting (root)
- CLAUDE.md Overview

## God Nodes (most connected - your core abstractions)
1. `LidarVNSync` - 77 edges
2. `timing_lidar` - 29 edges
3. `DopplerTurbulence` - 26 edges
4. `read_lidar` - 23 edges
5. `chunks` - 19 edges
6. `extract_sync_window()` - 16 edges
7. `refine_offset_20hz()` - 15 edges
8. `stare` - 14 edges
9. `NoaaDas` - 14 edges
10. `coarse_and_fine_lag()` - 13 edges

## Surprising Connections (you probably didn't know these)
- `Timing Offset Composition Equations` --references--> `chunk_lidar_datetimes()`  [EXTRACTED]
  .github/plans/OFFSET_EQUATIONS_AND_RESET_STRATEGY.md → lidar_vn_sync.jl
- `Sync Failure Modes & Prior Handling` --references--> `prior_from_history()`  [INFERRED]
  TODO.md → lidar_vn_sync.jl
- `Offset Composition Fix (single-residual rule)` --references--> `refine_offset_20hz()`  [EXTRACTED]
  .github/plans/EXECUTIVE_SUMMARY_OFFSET_FIX.md → lidar_vn_sync.jl
- `Timing Offset Composition Equations` --references--> `refine_offset_20hz()`  [EXTRACTED]
  .github/plans/OFFSET_EQUATIONS_AND_RESET_STRATEGY.md → lidar_vn_sync.jl
- `Graphify Workflow (Copilot Instructions)` --semantically_similar_to--> `graphify Mandatory Resource Rules (AGENTS.md)`  [INFERRED] [semantically similar]
  .github/copilot-instructions.md → AGENTS.md

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Four-Persona Agent Workflow (Orchestrator/Planner/Editor/Tester)** — github_agents_orchestrator_agent_orchestrator, github_agents_planner_agent_planner, github_agents_editor_agent_editor, github_agents_tester_agent_tester, agents_persona_routing [EXTRACTED 1.00]
- **Timing Offset Double-Count Bug: Diagnosis, Equations, Fix, and Failure-Mode Rationale** — github_plans_executive_summary_offset_fix_double_count_bug, github_plans_executive_summary_offset_fix_offset_composition_fix, github_plans_offset_equations_and_reset_strategy_offset_equations, github_plans_offset_equations_and_reset_strategy_prior_reset_strategy, todo_sync_failure_modes [INFERRED 0.85]
- **Structure-Function Dissipation Retrieval Pipeline** — readme_structure_function_fit, readme_epsilon_retrieval, production_add_qc_fit_stats_onepass, production_add_qc_se_epsilon_propagation, production_add_qc_epsilon_sentinels [INFERRED 0.85]

## Communities (51 total, 14 thin omitted)

### Community 0 - "read_lidar Core Chunking"
Cohesion: 0.08
Nodes (32): all_chunks(), all_gaps(), all_start_end_indices(), chunken(), chunkst(), compute_mdv_snr_mean(), get_all_file_start_end_idxs(), get_daily_meanuv() (+24 more)

### Community 1 - "LidarVNSync Core"
Cohesion: 0.07
Nodes (24): ensure_chunk_loaded!(), fcn3(), chunks, Dates, DSP, FFTW, Interpolations, JLD2 (+16 more)

### Community 2 - "Legacy Turbulence Cleanup (old/)"
Cohesion: 0.06
Nodes (18): chunks, Dates, DSP, FFTW, Interpolations, JLD2, NCDatasets, NoaaDas (+10 more)

### Community 3 - "NoaaDas Readers"
Cohesion: 0.16
Nodes (27): das_dict(), DasFps, DasGps, DasScs, declat(), declon(), get_das_filenames(), get_das_pathfiles() (+19 more)

### Community 4 - "Timing Offset Analysis"
Cohesion: 0.10
Nodes (21): ba(), cross_correlation(), dtregress(), DVM(), find_lags_iterative(), gps2dt(), gps2utc(), gpstime2gpsvndt() (+13 more)

### Community 5 - "DopplerTurbulence Module"
Cohesion: 0.08
Nodes (18): Module Consolidation into DopplerTurbulence.jl, DisplacementUncertainty, allcross(), displacements(), DopplerTurbulence, Dates, DSP, FFTW (+10 more)

### Community 6 - "ABL Working Group & Production QC"
Cohesion: 0.11
Nodes (21): Devmi, Entrainment and Cold Pools from Scalar Budgets, Jayesh, Lidar Backscatter BL Height (Haar Wavelet Gradient Detection), Mixing Depth from Lidar (Next Step), Simon, Epsilon Conversion epsilon(A), C2ll=2.0 (QC doc), Epsilon Sentinel Coding (QC doc) (+13 more)

### Community 7 - "Sync API Docs & Offset Sentinel"
Cohesion: 0.15
Nodes (18): LidarVNSync Core API Conventions (AGENTS.md), Sync-Quality Gates (AGENTS.md), LidarVNSync Core API Conventions (CLAUDE.md), Sync-Quality Gates (CLAUDE.md), Prior-Offset Double-Count Bug, Offset Composition Fix (single-residual rule), Prior Reset Strategy (carry_forward / reset_each_file / reset_on_bad_quality), append_nan_offset_log() (+10 more)

### Community 8 - "read_lidar Offset Fitting"
Cohesion: 0.16
Nodes (12): chunks, dt_to_chunkind(), findindices(), finelagcov(), fit_offset(), get_mean_uv_chunk!(), goodcov(), indavg() (+4 more)

### Community 9 - "Sync Window Extraction"
Cohesion: 0.20
Nodes (17): diag_array(), diagnostic_example_chunks(), diagnostic_single_window(), dt_seconds(), extract_sync_window(), fill_short_nan_gaps(), finite_overlap(), finite_overlap_corr() (+9 more)

### Community 10 - "VectorNav to JLD2 Conversion"
Cohesion: 0.13
Nodes (9): ShipPosmv, Dates, Interpolations, JLD2, NoaaDas, Pkg, read_vecnav, Revise (+1 more)

### Community 11 - "Chunk Timing Tests"
Cohesion: 0.17
Nodes (8): chunken(), chunkst(), Dates, read_lidar, Revise, nextchunki(), prevchunki(), thisj()

### Community 12 - "Hmix/Epsilon Computation"
Cohesion: 0.24
Nodes (12): compute_hmix_for_matrix(), extract_epsilon_and_meta(), hmix_first_crossing(), is_valid_eps(), Dates, JLD2, Pkg, Statistics (+4 more)

### Community 13 - "Deprecated Reader (2024-08-07)"
Cohesion: 0.21
Nodes (8): get_daily_meanuv(), Dates, Interpolations, JLD2, NCDatasets, read_lidar, read_streamlinexr_head(), read_streamlinexr_stare!()

### Community 14 - "Deprecated Reader (buggy offsets)"
Cohesion: 0.21
Nodes (8): get_daily_meanuv(), Dates, Interpolations, JLD2, NCDatasets, read_lidar, read_streamlinexr_head(), read_streamlinexr_stare!()

### Community 15 - "Lidar Index Building"
Cohesion: 0.20
Nodes (5): build_lidar_index(), Dates, JLD2, lidar_index, LidarIndex

### Community 16 - "Chunk Failure Scan"
Cohesion: 0.20
Nodes (7): Dates, FFTW, LidarVNSync, NCDatasets, Pkg, Printf, Statistics

### Community 17 - "Pitch/Roll Uncertainty & Cleanup Plan"
Cohesion: 0.24
Nodes (8): Import Organization (DisplacementUncertainty via DopplerTurbulence), Cleanup Plan After Pitch/Roll Uncertainty Implementation, Test Cell Reorganization (TEST1-3 to modules/other notebook), displacement_variance_pitch_roll(), DisplacementUncertainty, estimate_pitch_roll_uncertainty(), Statistics, propagate_rho_uncertainty()

### Community 18 - "Rain Despike Ideas (codeIdeas/foo.jl)"
Cohesion: 0.31
Nodes (7): estimate_heave_ts(), iterative_despike_rain(), Statistics, nanmedian(), rain_mask_vertical(), running_median(), subtract_rain_layers!()

### Community 19 - "FFT Lag/Xcorr & Offset Equations"
Cohesion: 0.31
Nodes (10): Timing Offset Composition Equations, analytic_envelope_fft(), coarse_and_fine_lag(), cosine_edge_mask(), detrend_center(), fft_bandpass(), fft_xcorr_lag(), iterative_coarse_lag() (+2 more)

### Community 20 - "POSMV/Gyro Readers"
Cohesion: 0.39
Nodes (9): Base.Iterators, get_nav_file(), get_posmv_file(), itr_expand(), read_gyro_data(), read_gyro_dict(), read_pashr_data(), read_pashr_dict() (+1 more)

### Community 21 - "Legacy Turbulence Binning/Plot"
Cohesion: 0.29
Nodes (8): binavg(), equal_bin(), f(), get_mdv(), indavg(), pcolor_lidar_stare(), read_stare_chunk(), remove_mdv()

### Community 22 - "Save Lidar Index Script"
Cohesion: 0.25
Nodes (6): JLD2, lidar_index, NCDatasets, Pkg, read_lidar, Revise

### Community 23 - "Epsilon JLD2 to NetCDF"
Cohesion: 0.29
Nodes (4): Dates, JLD2, NCDatasets, Pkg

### Community 24 - "FFT Derivative Test (old/)"
Cohesion: 0.33
Nodes (6): f_taper(), FFTW, Pkg, Revise, taper(), Plots

### Community 25 - "Legacy D2(rho) Structure Function"
Cohesion: 0.29
Nodes (7): allcross(), anom(), D2_rho_stare(), displacements(), lidarindices(), rng(), uniquepairs()

### Community 26 - "Save Lidar Datetime Script"
Cohesion: 0.29
Nodes (6): Dates, JLD2, NCDatasets, Pkg, read_lidar, Revise

### Community 27 - "Rain Flag Scan Idea"
Cohesion: 0.33
Nodes (4): FFTW, Statistics, ImageFiltering, LinearAlgebra

### Community 28 - "Legacy Offset Covariance Calc"
Cohesion: 0.33
Nodes (6): findindices(), good(), offset_cov(), offset_range_covs(), offset_subset(), sync_offset()

### Community 29 - "Save Hourly Stare NetCDF"
Cohesion: 0.33
Nodes (4): Dates, Pkg, Printf, read_lidar

### Community 30 - "Copilot 4-Persona Agents (legacy)"
Cohesion: 1.00
Nodes (5): Agent Init Table, Editor Agent, Orchestrator Agent, Planner Agent, Tester Agent

### Community 31 - "graphify Rules (Docs Cross-Ref)"
Cohesion: 0.50
Nodes (5): graphify Mandatory Resource Rules (AGENTS.md), AGENTS.md - Agent Routing & Workflow, graphify Rules (CLAUDE.md), Other Agent Contexts Note (legacy docs deprecated for Claude Code), Graphify Workflow (Copilot Instructions)

### Community 32 - "Iterative Heave/Rain Idea"
Cohesion: 0.60
Nodes (4): infer_heave_from_stare(), Statistics, local_vertical_stats(), weighted_median()

### Community 33 - "Lidar UV Script"
Cohesion: 0.40
Nodes (4): NCDatasets, Pkg, PythonPlot, Statistics

### Community 34 - "Chunk Alignment & Loading"
Cohesion: 0.67
Nodes (4): check_chunk_alignment_contract(), chunk_lidar_datetimes(), ensure_chunk_loaded_nc!(), load_chunk_for_dissipation()

### Community 35 - "Sync Context Setup"
Cohesion: 0.67
Nodes (4): init_periodic_beams(), load_lidar_indices_and_files(), setup_sync_context(), setup_sync_context_nc()

### Community 36 - "Motion Timing Doc (Leg1/Leg2)"
Cohesion: 0.50
Nodes (4): 18s GPS-UTC Leapsecond Offset (May-June 2024), Leg 1 VN-POSMV Clock Sync (13-hr precession), Leg 2 VN-POSMV Sync (0.6s lead, discretization artifact), Motion Timing Report (VectorNav/POSMV)

## Knowledge Gaps
- **144 isolated node(s):** `check-graphify-deps.sh script`, `run-in-lidar-env.sh script`, `bash`, `Revise`, `Pkg` (+139 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 248 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **14 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `LidarVNSync` connect `LidarVNSync Core` to `Chunk Alignment & Loading`, `Sync Context Setup`, `Sync API Docs & Offset Sentinel`, `Sync Window Extraction`, `FFT Lag/Xcorr & Offset Equations`?**
  _High betweenness centrality (0.024) - this node is a cross-community bridge._
- **Why does `Module Consolidation into DopplerTurbulence.jl` connect `DopplerTurbulence Module` to `Pitch/Roll Uncertainty & Cleanup Plan`, `ABL Working Group & Production QC`?**
  _High betweenness centrality (0.008) - this node is a cross-community bridge._
- **What connects `check-graphify-deps.sh script`, `run-in-lidar-env.sh script`, `bash` to the rest of the system?**
  _144 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `read_lidar Core Chunking` be split into smaller, more focused modules?**
  _Cohesion score 0.08170731707317073 - nodes in this community are weakly interconnected._
- **Should `LidarVNSync Core` be split into smaller, more focused modules?**
  _Cohesion score 0.06554621848739496 - nodes in this community are weakly interconnected._
- **Should `Legacy Turbulence Cleanup (old/)` be split into smaller, more focused modules?**
  _Cohesion score 0.058823529411764705 - nodes in this community are weakly interconnected._
- **Should `Timing Offset Analysis` be split into smaller, more focused modules?**
  _Cohesion score 0.10344827586206896 - nodes in this community are weakly interconnected._