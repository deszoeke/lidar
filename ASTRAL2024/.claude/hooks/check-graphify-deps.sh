#!/bin/bash
# SessionStart hook: warn if the lidar micromamba env is missing packages
# the graphify MCP server needs (mcp, graphify). This repo is synced via
# git across multiple computers, so envs can drift.

missing=""
for pkg in mcp graphify; do
  bash "$(dirname "$0")/run-in-lidar-env.sh" python -c "import $pkg" >/dev/null 2>&1 || missing="$missing $pkg"
done

if [ -n "$missing" ]; then
  context="graphify: the 'lidar' micromamba env on this machine is missing python package(s):$missing. Tell the user and offer to run: micromamba run -n lidar pip install$missing (or micromamba install -n lidar -c conda-forge$missing). The graphify MCP server will not start without them; file-based graphify (graphify-out/GRAPH_REPORT.md, graphify update .) still works in the meantime."
  python3 -c "import json,sys; print(json.dumps({'hookSpecificOutput':{'hookEventName':'SessionStart','additionalContext':sys.argv[1]}}))" "$context"
fi
