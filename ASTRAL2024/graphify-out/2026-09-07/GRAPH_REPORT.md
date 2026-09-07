# Graph Report - ASTRAL2024  (2026-09-07)

## Corpus Check
- 43 files · ~44,293 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 681 nodes · 1026 edges · 44 communities (40 shown, 4 thin omitted)
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `eca245c9`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- LidarVNSync
- JLD2
- ASTRAL2024/old/lidar_turbulence_cleanup.jl
- read_lidar
- NoaaDas
- timing_lidar
- DopplerTurbulence
- Statistics
- chunks
- time_tests.jl
- old/lidar_turbulence_cleanup.jl
- lidar_index
- compute_hmix_epsilon.jl
- Timing Offset Equations & Prior Reset Strategy
- Agent Routing & Workflow
- vectornav2jld.jl
- Executive Summary: Timing Offset Double-Count Bug Fix
- Cleanup Plan After Pitch/Roll Uncertainty Implementation
- chunk_failure_scan.jl
- read_lidar
- read_lidar
- codeIdeas/foo.jl
- Dates
- read_lidar
- read_lidar
- old/fft_derivative_test.jl
- Turbulence processing
- save_lidar_dt.jl
- save_lidar_index.jl
- epsi_jld2nc.jl
- Editor
- Orchestrator
- Planner
- Tester
- Leg 1
- Production Run Additions
- codeIdeas/iterative_heave_rain.jl
- lidar_uv.jl
- ABLwh20260409.md
- save_stare_mdv_sync_pass2.jl
- CLAUDE.md
- extend_pycall_missing.jl
- copilot-instructions.md

## God Nodes (most connected - your core abstractions)
1. `LidarVNSync` - 91 edges
2. `DopplerTurbulence` - 34 edges
3. `timing_lidar` - 34 edges
4. `read_lidar` - 28 edges
5. `Dates` - 23 edges
6. `chunks` - 22 edges
7. `NoaaDas` - 17 edges
8. `extract_sync_window()` - 16 edges
9. `stare` - 16 edges
10. `lidar_index` - 13 edges

## Surprising Connections (you probably didn't know these)
- None detected - all connections are within the same source files.

## Import Cycles
- None detected.

## Communities (44 total, 4 thin omitted)

### Community 0 - "LidarVNSync"
Cohesion: 0.05
Nodes (68): analytic_envelope_fft(), append_nan_offset_log(), backward_jump_robustness(), check_chunk_alignment_contract(), chunk_lidar_datetimes(), coarse_and_fine_lag(), cosine_edge_mask(), detrend_center() (+60 more)

### Community 1 - "JLD2"
Cohesion: 0.24
Nodes (8): to_cf(), f_taper(), taper(), JLD2, NCDatasets, Pkg, Revise, load_file_beam_inds()

### Community 2 - "ASTRAL2024/old/lidar_turbulence_cleanup.jl"
Cohesion: 0.07
Nodes (27): cat_dicts, allcross(), anom(), binavg(), D2_rho_stare(), displacements(), equal_bin(), f() (+19 more)

### Community 3 - "read_lidar"
Cohesion: 0.08
Nodes (32): all_chunks(), all_gaps(), all_start_end_indices(), chunken(), chunkst(), compute_mdv_snr_mean(), get_all_file_start_end_idxs(), get_daily_meanuv() (+24 more)

### Community 4 - "NoaaDas"
Cohesion: 0.11
Nodes (38): Base.Iterators, das_dict, flatten, das_dict(), DasFps, DasGps, DasScs, declat() (+30 more)

### Community 5 - "timing_lidar"
Cohesion: 0.09
Nodes (22): StatsBase, ba(), cross_correlation(), dtregress(), DVM(), find_lags_iterative(), gps2dt(), gps2utc() (+14 more)

### Community 6 - "DopplerTurbulence"
Cohesion: 0.08
Nodes (17): DisplacementUncertainty, allcross(), displacements(), DopplerTurbulence, Dates, DSP, FFTW, Interpolations (+9 more)

### Community 7 - "Statistics"
Cohesion: 0.08
Nodes (18): FFTW, Statistics, DisplacementUncertainty, Statistics, FFTW, estimate_heave_ts(), iterative_despike_rain(), nanmedian() (+10 more)

### Community 8 - "chunks"
Cohesion: 0.16
Nodes (12): chunks, dt_to_chunkind(), findindices(), finelagcov(), fit_offset(), get_mean_uv_chunk!(), goodcov(), indavg() (+4 more)

### Community 9 - "time_tests.jl"
Cohesion: 0.24
Nodes (14): chunken(), chunkst(), get_end_fileidx(), get_stare_files(), get_start_fileidx(), isend(), isstart(), Dates (+6 more)

### Community 10 - "old/lidar_turbulence_cleanup.jl"
Cohesion: 0.05
Nodes (40): allcross(), anom(), binavg(), D2_rho_stare(), displacements(), equal_bin(), f(), findindices() (+32 more)

### Community 11 - "lidar_index"
Cohesion: 0.18
Nodes (5): build_lidar_index(), Dates, JLD2, lidar_index, LidarIndex

### Community 12 - "compute_hmix_epsilon.jl"
Cohesion: 0.30
Nodes (12): compute_hmix_for_matrix(), extract_epsilon_and_meta(), hmix_first_crossing(), is_valid_eps(), Dates, JLD2, Pkg, Statistics (+4 more)

### Community 13 - "Timing Offset Equations & Prior Reset Strategy"
Cohesion: 0.10
Nodes (20): Change 1: `run_sequential_offsets()` signature, Change 2: `prior_from_history()` logic, Change 3: Track file index in `run_sequential_offsets()`, Change 4: Fix the double-count bug in `process_sync_data()`, Correct Equations (What We Want), Current Implementation (Broken), Equation 1: Lidar time synchronization, Equation 2: 1 Hz coarse + fine lag stage (+12 more)

### Community 14 - "Agent Routing & Workflow"
Cohesion: 0.12
Nodes (16): Agent Descriptions & When LLMs Must Use Them, Agent Init, Agent Routing & Workflow, Caveman Mode (Active), Editor, Environment, graphify (Mandatory Resource), Key data structures (+8 more)

### Community 15 - "vectornav2jld.jl"
Cohesion: 0.15
Nodes (15): read_vecnav_dict, ShipPosmv, anom(), cat_dicts(), f(), Dates, Interpolations, JLD2 (+7 more)

### Community 16 - "Executive Summary: Timing Offset Double-Count Bug Fix"
Cohesion: 0.13
Nodes (14): After (New Code):, Before (Old Code):, Code Changes, Configurable Prior Reset (if cascade failures are a problem), Executive Summary: Timing Offset Double-Count Bug Fix, Files Updated, How to Verify, Next Steps (Optional, Not Critical) (+6 more)

### Community 17 - "Cleanup Plan After Pitch/Roll Uncertainty Implementation"
Cohesion: 0.14
Nodes (13): Additional Issues Found During Testing, Cleanup Plan After Pitch/Roll Uncertainty Implementation, Implementation Order, Import Organization, Issue: Check other dependents of DopplerTurbulence.jl, Issue: DisplacementUncertainty imported mid-notebook, Issue: Methods duplicated between DopplerTurbulence.jl and notebook, Issue: TEST 1 uses VN-aligned data, doesn't belong in displacement_uncertainty (+5 more)

### Community 18 - "chunk_failure_scan.jl"
Cohesion: 0.18
Nodes (10): band_coherence_welch_5_20s(), Dates, FFTW, LidarVNSync, NCDatasets, Pkg, Printf, Statistics (+2 more)

### Community 19 - "read_lidar"
Cohesion: 0.21
Nodes (8): get_daily_meanuv(), Dates, Interpolations, JLD2, NCDatasets, read_lidar, read_streamlinexr_head(), read_streamlinexr_stare!()

### Community 20 - "read_lidar"
Cohesion: 0.21
Nodes (8): get_daily_meanuv(), Dates, Interpolations, JLD2, NCDatasets, read_lidar, read_streamlinexr_head(), read_streamlinexr_stare!()

### Community 21 - "codeIdeas/foo.jl"
Cohesion: 0.31
Nodes (7): estimate_heave_ts(), iterative_despike_rain(), Statistics, nanmedian(), rain_mask_vertical(), running_median(), subtract_rain_layers!()

### Community 22 - "Dates"
Cohesion: 0.22
Nodes (7): Dates, read_streamlinexr_stare, find_hpl_files(), Dates, Pkg, Printf, read_lidar

### Community 23 - "read_lidar"
Cohesion: 0.31
Nodes (5): Interpolations, get_daily_meanuv(), read_lidar, read_streamlinexr_head(), read_streamlinexr_stare!()

### Community 24 - "read_lidar"
Cohesion: 0.36
Nodes (4): get_daily_meanuv(), read_lidar, read_streamlinexr_head(), read_streamlinexr_stare!()

### Community 25 - "old/fft_derivative_test.jl"
Cohesion: 0.33
Nodes (6): f_taper(), FFTW, Pkg, Revise, taper(), Plots

### Community 26 - "Turbulence processing"
Cohesion: 0.29
Nodes (6): Data to read, Fit statistics and postprocessing, Production workflow, Routines for reading and processing Doppler wind lidar, Saved production outputs, Turbulence processing

### Community 27 - "save_lidar_dt.jl"
Cohesion: 0.29
Nodes (6): Dates, JLD2, NCDatasets, Pkg, read_lidar, Revise

### Community 28 - "save_lidar_index.jl"
Cohesion: 0.29
Nodes (6): JLD2, lidar_index, NCDatasets, Pkg, read_lidar, Revise

### Community 29 - "epsi_jld2nc.jl"
Cohesion: 0.33
Nodes (5): Dates, JLD2, NCDatasets, Pkg, to_f64()

### Community 30 - "Editor"
Cohesion: 0.33
Nodes (5): Code Style Context, Core Responsibilities, Editor, Output Format (Caveman Clause), When to Use This Agent

### Community 31 - "Orchestrator"
Cohesion: 0.33
Nodes (5): Core Responsibilities, Implementation Style, Orchestrator, Output Format (Caveman Clause), When to Use This Agent

### Community 32 - "Planner"
Cohesion: 0.33
Nodes (5): Code Style Context, Core Responsibilities, Output Format (Caveman Clause), Planner, When to Use This Agent

### Community 33 - "Tester"
Cohesion: 0.33
Nodes (5): Best Practices, Core Responsibilities, Output Format (Caveman Clause), Tester, When to Use This Agent

### Community 34 - "Leg 1"
Cohesion: 0.33
Nodes (5): A Report gleaned from experiments in [timing_lidar.jl](timing_lidar.jl) and [lidar_turbulence.jl](lidar_turbulence.jl), Caption: (+) Positive lag means Vn start index is shifted forward, i.e., Vn signal lags POSMV. (Data in this shifted index is (moved backward and) aligned with the unshifted POSMV start.) (-) Negative lag means negative index of Vn lines up with start of POSMV epoch, i.e., Vn signal leads POSMV., Leg 1, Leg 2, part 1, Timing of motion data for the Halo Photonics lidar from VectorNav and POSMV

### Community 35 - "Production Run Additions"
Cohesion: 0.33
Nodes (5): Efficiency changes, Output variables saved, Postproduction utility, Principles of calculation, Production Run Additions

### Community 36 - "codeIdeas/iterative_heave_rain.jl"
Cohesion: 0.60
Nodes (4): infer_heave_from_stare(), Statistics, local_vertical_stats(), weighted_median()

### Community 37 - "lidar_uv.jl"
Cohesion: 0.40
Nodes (4): NCDatasets, Pkg, PythonPlot, Statistics

### Community 38 - "ABLwh20260409.md"
Cohesion: 0.50
Nodes (3): Agenda, EKAMSAT ABL Working Group, Lidar turbulence analysis - Simon

## Knowledge Gaps
- **199 isolated node(s):** `Revise`, `Pkg`, `Dates`, `Statistics`, `Interpolations` (+194 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **4 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Dates` connect `Dates` to `LidarVNSync`, `JLD2`, `ASTRAL2024/old/lidar_turbulence_cleanup.jl`, `read_lidar`, `NoaaDas`, `timing_lidar`, `DopplerTurbulence`, `chunks`, `time_tests.jl`, `lidar_index`, `compute_hmix_epsilon.jl`, `vectornav2jld.jl`, `chunk_failure_scan.jl`, `read_lidar`, `read_lidar`?**
  _High betweenness centrality (0.221) - this node is a cross-community bridge._
- **Why does `LidarVNSync` connect `LidarVNSync` to `JLD2`, `ASTRAL2024/old/lidar_turbulence_cleanup.jl`, `timing_lidar`, `DopplerTurbulence`, `Statistics`, `vectornav2jld.jl`, `chunk_failure_scan.jl`, `Dates`, `read_lidar`?**
  _High betweenness centrality (0.168) - this node is a cross-community bridge._
- **Why does `PyPlot` connect `ASTRAL2024/old/lidar_turbulence_cleanup.jl` to `old/lidar_turbulence_cleanup.jl`?**
  _High betweenness centrality (0.109) - this node is a cross-community bridge._
- **What connects `Revise`, `Pkg`, `Dates` to the rest of the system?**
  _199 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `LidarVNSync` be split into smaller, more focused modules?**
  _Cohesion score 0.05063291139240506 - nodes in this community are weakly interconnected._
- **Should `ASTRAL2024/old/lidar_turbulence_cleanup.jl` be split into smaller, more focused modules?**
  _Cohesion score 0.06871035940803383 - nodes in this community are weakly interconnected._
- **Should `read_lidar` be split into smaller, more focused modules?**
  _Cohesion score 0.08130081300813008 - nodes in this community are weakly interconnected._