# Agent Routing & Workflow

## When to Invoke Which Agent

**User asks for help?** → Invoke **Orchestrator** (always start here)
- Orchestrator decides the sequence and delegates to planner/editor/tester
- Orchestrator coordinates the workflow and reports findings to user

**User asks for a plan/architecture review/code analysis?** → Invoke **Planner** directly
- Use when: analyzing code, detecting bugs, assessing risks, reviewing impact
- Use when: designing implementation steps before major edits
- Planner reads graphify and proposes changes; hands off to editor for execution

**User asks to implement a plan or make a code edit?** → Invoke **Editor** 
- Use when: executing approved changes (from planner), minimal surgical edits
- Use when: notebook cell modifications, function rewrites, file creation
- Editor makes precise changes, preserves existing style

**User asks to test, validate, or diagnose?** → Invoke **Tester**
- Use when: writing short validation tests, smoke tests, diagnostics  
- Use when: checking schemas, key outputs, or reproducing failures
- Tester creates repeatable test code and runs short tests; user runs long jobs

---

## Agent Descriptions & When LLMs Must Use Them

### Orchestrator 
- **When to invoke:** User asks "help me with X" or provides general request
- **Role:** Routes to planner → editor → tester in correct sequence
- **Decides:** What order to work in, when planning vs. direct edit is best
- **Delegates:** Planner for analysis, Editor for changes, Tester for validation
- **OUTPUT:** Status updates, next-step recommendations, key findings summarized (not raw analysis)

### Planner
- **When to invoke:** User asks to "plan," "analyze," "review," "architect," or "check for bugs"
- **Role:** Creates implementation plans, detects regressions, identifies risks
- **Reads:** graphify-out/GRAPH_REPORT.md and graphify-out/graph.json for dependencies
- **Produces:** Phased steps, acceptance criteria, impact analysis, mitigation strategies
- **OUTPUT:** Clear plan + risk assessment + rollback points (numbered steps, not prose)

### Editor
- **When to invoke:** Orchestrator/Planner hands off a plan to execute
- **Role:** Makes minimal, targeted code changes; preserves style and APIs
- **Executes:** Only approved changes from Planner; no scope creep
- **Validates:** Quick smoke tests after edits
- **OUTPUT:** Confirmation of what changed + link to modified files (short, surgical summary)

### Tester
- **When to invoke:** Need validation, diagnostics, or test setup
- **Role:** Writes short repeatable tests; prepares long job commands for user
- **Creates:** Notebook cells, scripts, or terminal commands user can rerun
- **Does NOT:** Run long production jobs without explicit request
- **OUTPUT:** Test code + expected results + clear pass/fail criteria (executable artifacts)

---

## graphify (Mandatory Resource)

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

---

## Environment

This project runs under the `lidar` mamba environment. Before running any terminal commands, activate it:

```bash
source ~/miniforge3/etc/profile.d/mamba.sh && mamba activate lidar
```
