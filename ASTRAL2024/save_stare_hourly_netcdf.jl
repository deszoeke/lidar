"""
save_stare_hourly_netcdf.jl

Convert all Stare_116_YYYYMMDD_HH.hpl files under data/YYYYMMDD/ to
NetCDF-4 files in data/netcdf_stare/.

Usage (from the project root):
    julia save_stare_hourly_netcdf.jl              # convert all, skip existing
    julia save_stare_hourly_netcdf.jl --overwrite  # reconvert all
    julia save_stare_hourly_netcdf.jl 20240501     # one day only
"""

cd(@__DIR__)
using Pkg; Pkg.activate(".")

include("read_lidar.jl")
using .read_lidar: read_streamlinexr_stare,
                   compute_mdv_snr_mean, write_stare_hourly_netcdf

using Dates
using Printf

# ── configuration ──────────────────────────────────────────────────────────
const LIDAR_DATA_DIR = joinpath(@__DIR__, "data")
const NETCDF_OUT_DIR = joinpath(@__DIR__, "data", "netcdf_stare")
const SNR_THRESHOLD  = 1.03          # intensity > 1.03 ↔ SNR > 0.03

# ── command-line arguments ─────────────────────────────────────────────────
overwrite   = "--overwrite" in ARGS
day_filter  = filter(a -> !startswith(a, "-"), ARGS)  # bare date args

mkpath(NETCDF_OUT_DIR)

# ── discover HPL files ────────────────────────────────────────────────────
"""
Return all Stare_116_YYYYMMDD_HH.hpl paths under LIDAR_DATA_DIR,
optionally filtered to the requested yyyymmdd strings.
"""
function find_hpl_files(datadir::AbstractString;
                         day_filter::Vector{String} = String[])
    hpl_files = String[]
    for daydir in readdir(datadir; join=true)
        isdir(daydir) || continue
        dayname = basename(daydir)
        # skip if day_filter is given and this day is not listed
        if !isempty(day_filter) && !(dayname in day_filter)
            continue
        end
        for fname in readdir(daydir; join=true)
            if occursin(r"^Stare_116_\d{8}_\d{2}\.hpl$", basename(fname))
                push!(hpl_files, fname)
            end
        end
    end
    sort!(hpl_files)
    return hpl_files
end

hpl_files = find_hpl_files(LIDAR_DATA_DIR; day_filter)

n_total  = length(hpl_files)
n_done   = 0
n_skip   = 0
n_err    = 0

@printf("Found %d HPL files.  Output dir: %s\n", n_total, NETCDF_OUT_DIR)
overwrite && println("Mode: overwrite existing files.")

for (i, hpl_path) in enumerate(hpl_files)
    stem    = splitext(basename(hpl_path))[1]   # "Stare_116_YYYYMMDD_HH"
    nc_path = joinpath(NETCDF_OUT_DIR, stem * ".nc")

    if !overwrite && isfile(nc_path)
        n_skip += 1
        continue
    end

    try
        # ── read header and data in one call ─────────────────────────
        beams, header = read_streamlinexr_stare(hpl_path)
        header[:source_hpl] = hpl_path          # store provenance

        # ── write NetCDF ──────────────────────────────────────────────
        write_stare_hourly_netcdf(nc_path, beams, header;
                                  snr_threshold = SNR_THRESHOLD)
        n_done += 1

        if mod(i, 50) == 0 || i == n_total
            @printf("[%4d/%4d] done=%d skip=%d err=%d\n",
                    i, n_total, n_done, n_skip, n_err)
        end

    catch e
        n_err += 1
        @printf("ERROR  [%d/%d] %s\n  %s\n", i, n_total, hpl_path, e)
    end
end

@printf("\nFinished.  Converted=%d  Skipped=%d  Errors=%d\n",
        n_done, n_skip, n_err)
