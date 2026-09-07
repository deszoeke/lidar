## Style: Caveman Mode (default)

Use short, direct sentences. Prefer simple words over formal wording.
Keep status updates compact and action-first: status, key result, next step.
Avoid long prose and unnecessary theory. If the user asks for deeper detail, give it without
dropping Caveman Mode for later turns.

## Workflow

No persona relay (no Orchestrator/Planner/Editor/Tester hand-off). Use Claude Code's native
loop instead:
- **Planning** — use plan mode for nontrivial or risky edits. Plan mode does not know to use
  graphify on its own: before proposing a plan, or answering any architecture/impact question,
  read `graphify-out/GRAPH_REPORT.md` (god nodes, communities, blast radius) instead of grepping
  raw files first.
- **Editing** — done directly in the main session. Minimal, surgical changes; preserve existing
  APIs and style; no unrelated refactors or formatting churn.
- **Testing** — write short, repeatable validation code (notebook cells or scripts) the user can
  rerun and edit. Do not launch long production runs without explicit request; break long runs
  into short testable pieces first.
- **Context hygiene** — if a task is about to involve noisy, repetitive exploration (reading many
  files, large grep dumps, long log tails), say so and suggest delegating it to a subagent
  (e.g. Explore, or a forked/general-purpose agent) instead of doing it inline and bloating this
  session's context.

## graphify

This project has a graphify knowledge graph at `graphify-out/`, queryable two ways:
- **MCP** (`.mcp.json`, server `graphify`) — live queries against the graph.
- **Static files** — `graphify-out/GRAPH_REPORT.md` (god nodes, communities, gaps) and
  `graphify-out/graph.json`.

Rules:
- Before answering architecture or codebase questions, read `graphify-out/GRAPH_REPORT.md` for
  god nodes and community structure.
- Before proposing any nontrivial edit plan, consult the graph (MCP or `GRAPH_REPORT.md`/
  `graph.json`) to identify dependencies, central modules, and transitive impact.
- If `graphify-out/wiki/index.md` exists, navigate it instead of reading raw files, unless a raw
  read is needed to confirm implementation details.
- If MCP graph tools appear stale but `graphify-out/graph.json` or `GRAPH_REPORT.md` is fresh,
  trust the on-disk artifacts and continue.
- After modifying code files in this session, run `graphify update .` to keep the graph current
  (AST-only, no API cost).
- After major structural edits, re-read `GRAPH_REPORT.md` before summarizing architecture or
  dependency impact.

## Environment

This project runs under the `lidar` micromamba/conda env. Activate before running terminal
commands:

```bash
# micromamba (preferred)
eval "$(micromamba shell hook --shell zsh)" && micromamba activate lidar

# mamba / conda fallback
source ~/miniforge3/etc/profile.d/mamba.sh && mamba activate lidar
```

Use whichever succeeds; verify the prompt/env shows `lidar` before proceeding.

## LidarVNSync Core API Conventions

### Key data structures

| Symbol | Type | Meaning |
|--------|------|---------|
| `icvn` | `UnitRange{Int}` | Chunk indices covered by VectorNav data (from `setup_sync_context`) |
| `Env.ists[ic]` | `Int` | Index of the first lidar beam in chunk `ic` |
| `Env.dtime[Env.ists[ic]]` | `DateTime` | Wall-clock start time of chunk `ic` |
| `Env.dtime[Env.iens[ic]]` | `DateTime` | Wall-clock end time of chunk `ic` |

### Selecting chunks by date — correct idiom

```julia
# correct: lambda contains the full predicate, collect(icvn) is the collection
ic_day1 = filter(ic -> Date(Env.dtime[Env.ists[ic]]) == Date(2024, 4, 28), collect(icvn))
```

Common mistake — a stray `) == Date(...)` outside the lambda produces a parse error:
```julia
# wrong (syntax error):
ic_day1 = filter(ic -> Date(Env.dtime[Env.ists[ic]]) == Date(2024,4,28)) == Date(2024,4,28), collect(icvn))
```

### Key function signatures

```julia
# Run sync over an explicit chunk list (ic_list=nothing → all icvn)
LidarVNSync.write_daily_mdv_vn2!(;
    out_dir  = "data/vn_sync_chunk_daily",
    ic_list  = ic_day1,   # Vector{Int} or nothing
    overwrite = false,
)

# Sequential offset computation (returns NamedTuple result)
result = LidarVNSync.process_sync_data(beams, Env, Vn, UV, ic_list; ntop=nz)

# Single-chunk window extraction (stateful: pass the same `state` across calls)
win = LidarVNSync.extract_sync_window(beams, Env, state, Vn, UV, ic; ntop=nz, nc_dir="./data/netcdf_stare")
```

### Sync-quality gates

A chunk is skipped (offset → sentinel `−9999`) when:
- `win.vn_coverage < 0.5` — VN 20 Hz data covers < 50% of the window
- `win.vn2_xcorr_nan_frac > 0.15` — too many NaN samples in the interpolated VN signal

Downstream code must guard: `isfinite(offset_s) && offset_s != LidarVNSync.OFFSET_SENTINEL_S`.

## Other agent contexts in this repo

`AGENTS.md`, `.github/agents/`, `.github/copilot-instructions.md`, and `.codex/hooks.json` are
kept as a fallback/cross-LLM safety net (Copilot/Codex), not actively maintained for Claude Code.
Don't assume Claude Code follows their persona-routing content — this file is authoritative for
Claude Code.
