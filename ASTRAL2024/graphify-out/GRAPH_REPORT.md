# Graph Report - /home/deszoeks/projects/lidar/ASTRAL2024  (2026-09-04)

## Corpus Check
- 25 files · ~66,453 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 350 nodes · 641 edges · 14 communities detected
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
- [[_COMMUNITY_Community 10|Community 10]]
- [[_COMMUNITY_Community 11|Community 11]]
- [[_COMMUNITY_Community 12|Community 12]]
- [[_COMMUNITY_Community 13|Community 13]]

## God Nodes (most connected - your core abstractions)
1. `LidarVNSync` - 75 edges
2. `timing_lidar` - 29 edges
3. `DopplerTurbulence` - 25 edges
4. `Dates` - 23 edges
5. `read_lidar` - 23 edges
6. `chunks` - 19 edges
7. `extract_sync_window()` - 16 edges
8. `NoaaDas` - 14 edges
9. `stare` - 14 edges
10. `refine_offset_20hz()` - 13 edges

## Surprising Connections (you probably didn't know these)
- None detected - all connections are within the same source files.

## Communities

### Community 0 - "Community 0"
Cohesion: 0.07
Nodes (53): analytic_envelope_fft(), append_nan_offset_log(), backward_jump_robustness(), check_chunk_alignment_contract(), chunk_lidar_datetimes(), coarse_and_fine_lag(), cosine_edge_mask(), detrend_center() (+45 more)

### Community 1 - "Community 1"
Cohesion: 0.06
Nodes (24): Dates, f_taper(), taper(), Interpolations, JLD2, NCDatasets, Pkg, Plots (+16 more)

### Community 2 - "Community 2"
Cohesion: 0.07
Nodes (27): cat_dicts, allcross(), anom(), binavg(), D2_rho_stare(), displacements(), equal_bin(), f() (+19 more)

### Community 3 - "Community 3"
Cohesion: 0.11
Nodes (25): all_chunks(), all_gaps(), all_start_end_indices(), chunken(), chunkst(), compute_mdv_snr_mean(), get_all_file_start_end_idxs(), get_daily_meanuv() (+17 more)

### Community 4 - "Community 4"
Cohesion: 0.16
Nodes (25): das_dict, das_dict(), DasFps, DasGps, DasScs, declat(), declon(), get_das_filenames() (+17 more)

### Community 5 - "Community 5"
Cohesion: 0.12
Nodes (18): StatsBase, ba(), cross_correlation(), dtregress(), DVM(), find_lags_iterative(), gps2dt(), gps2utc() (+10 more)

### Community 6 - "Community 6"
Cohesion: 0.12
Nodes (9): allcross(), displacements(), DopplerTurbulence, lidarindices(), rng(), trigs(), uniquepairs(), wtrue() (+1 more)

### Community 7 - "Community 7"
Cohesion: 0.12
Nodes (8): DisplacementUncertainty, FFTW, ImageFiltering, infer_heave_from_stare(), local_vertical_stats(), weighted_median(), LinearAlgebra, Statistics

### Community 8 - "Community 8"
Cohesion: 0.18
Nodes (11): chunks, dt_to_chunkind(), findindices(), finelagcov(), fit_offset(), get_mean_uv_chunk!(), goodcov(), indavg() (+3 more)

### Community 9 - "Community 9"
Cohesion: 0.23
Nodes (5): chunken(), chunkst(), nextchunki(), prevchunki(), thisj()

### Community 10 - "Community 10"
Cohesion: 0.39
Nodes (9): flatten, get_nav_file(), get_posmv_file(), itr_expand(), read_gyro_data(), read_gyro_dict(), read_pashr_data(), read_pashr_dict() (+1 more)

### Community 11 - "Community 11"
Cohesion: 0.25
Nodes (3): build_lidar_index(), lidar_index, LidarIndex

### Community 12 - "Community 12"
Cohesion: 0.42
Nodes (8): compute_hmix_for_matrix(), extract_epsilon_and_meta(), hmix_first_crossing(), is_valid_eps(), list_epsilon_files(), main(), process_file(), rng_height()

### Community 13 - "Community 13"
Cohesion: 0.36
Nodes (6): estimate_heave_ts(), iterative_despike_rain(), nanmedian(), rain_mask_vertical(), running_median(), subtract_rain_layers!()

## Knowledge Gaps
- **9 isolated node(s):** `das_dict`, `flatten`, `read_stare_chunk`, `read_streamlinexr_stare`, `LinearAlgebra` (+4 more)
  These have ≤1 connection - possible missing edges or undocumented components.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Dates` connect `Community 1` to `Community 0`, `Community 2`, `Community 3`, `Community 4`, `Community 5`, `Community 6`, `Community 7`, `Community 8`, `Community 9`, `Community 10`, `Community 11`, `Community 12`?**
  _High betweenness centrality (0.461) - this node is a cross-community bridge._
- **Why does `LidarVNSync` connect `Community 0` to `Community 1`, `Community 2`, `Community 5`, `Community 6`, `Community 7`?**
  _High betweenness centrality (0.349) - this node is a cross-community bridge._
- **Why does `Statistics` connect `Community 7` to `Community 0`, `Community 1`, `Community 2`, `Community 5`, `Community 6`, `Community 8`, `Community 12`, `Community 13`?**
  _High betweenness centrality (0.185) - this node is a cross-community bridge._
- **What connects `das_dict`, `flatten`, `read_stare_chunk` to the rest of the system?**
  _9 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.07 - nodes in this community are weakly interconnected._
- **Should `Community 1` be split into smaller, more focused modules?**
  _Cohesion score 0.06 - nodes in this community are weakly interconnected._
- **Should `Community 2` be split into smaller, more focused modules?**
  _Cohesion score 0.07 - nodes in this community are weakly interconnected._