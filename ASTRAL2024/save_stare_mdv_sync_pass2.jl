"""
save_stare_mdv_sync_pass2.jl

Second pass over hourly stare NetCDF files that appends `mdv_sync(time)`.
`mdv_sync` is computed with the same chunk sync path used by
`extract_sync_window`/`valid_chunk_mdv`, so it matches sync-loop mdv.

Usage:
    julia save_stare_mdv_sync_pass2.jl
    julia save_stare_mdv_sync_pass2.jl --overwrite
    julia save_stare_mdv_sync_pass2.jl --ntop=80
"""

cd(@__DIR__)
using Pkg; Pkg.activate(".")

include("lidar_vn_sync.jl")
using .LidarVNSync

overwrite = "--overwrite" in ARGS
ntop = 80
for a in ARGS
    startswith(a, "--ntop=") || continue
    global ntop = parse(Int, split(a, "=", limit=2)[2])
end

res = LidarVNSync.write_mdv_sync_pass2!(;
    nc_dir=joinpath("data", "netcdf_stare"),
    ntop=ntop,
    overwrite=overwrite,
    log_every=200,
)

println("Pass-2 complete: ", res)
