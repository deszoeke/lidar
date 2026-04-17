## graphify

This project has a graphify knowledge graph at graphify-out/.

Rules:
- Treat graphify as a core skill for this repository, not an optional tool.
- Before answering architecture or codebase questions, read graphify-out/GRAPH_REPORT.md for god nodes, communities, and likely blast radius.
- Before proposing any nontrivial edit plan, consult graphify-out/GRAPH_REPORT.md and graphify-out/graph.json to identify dependencies, central modules, and transitive impact.
- If graphify-out/wiki/index.md exists, navigate it instead of reading raw files unless a raw-file read is needed to confirm implementation details.
- Prefer graph-informed planning: use the graph to identify directly affected modules, likely downstream consumers, and focused validation targets.
- If MCP graph tools appear stale but graphify-out/graph.json or graphify-out/GRAPH_REPORT.md is fresh, trust the on-disk graph artifacts and continue.
- After modifying code files in this session, run `graphify update .` to keep the graph current (AST-only, no API cost).
- After major structural edits, re-read graphify-out/GRAPH_REPORT.md before summarizing architecture or dependency impact.
