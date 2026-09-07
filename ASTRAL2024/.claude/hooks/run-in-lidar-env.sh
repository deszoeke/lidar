#!/bin/bash
# Run a command in the `lidar` env, preferring micromamba, falling back to
# mamba/conda (miniforge3 layout) if micromamba isn't available on this
# machine. This repo is synced via git across multiple computers whose
# mamba setups can differ.

if command -v micromamba >/dev/null 2>&1; then
  exec micromamba run -n lidar "$@"
elif [ -f ~/miniforge3/etc/profile.d/conda.sh ]; then
  source ~/miniforge3/etc/profile.d/conda.sh
  exec mamba run -n lidar "$@"
else
  echo "run-in-lidar-env.sh: neither micromamba nor mamba/conda (miniforge3) found" >&2
  exit 1
fi
