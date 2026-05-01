using Revise
using Pkg; Pkg.activate(".")

using JLD2
using NCDatasets

include("./lidar_index.jl")
using .lidar_index

include("./read_lidar.jl")
using .read_lidar

function load_file_beam_inds(files)
    if isfile("file_beam_inds.jld2")
        fileinds = load("file_beam_inds.jld2")
        file_ibeam_start = Int.(fileinds["bigind_file_start"])
        file_ibeam_end = Int.(fileinds["bigind_file_end"])
        if length(file_ibeam_start) == length(files)
            return file_ibeam_start, file_ibeam_end
        end
        @warn "file_beam_inds.jld2 has $(length(file_ibeam_start)) files but directory has $(length(files)) files; recomputing."
    end

    nfiles = length(files)
    nbeams = zeros(Int, nfiles)

    for (i, file) in enumerate(files)
        NCDataset(file, "r") do ds
            nbeams[i] = length(ds["time"])
        end
    end

    file_ibeam_end = cumsum(nbeams)
    file_ibeam_start = [1; file_ibeam_end[1:end-1] .+ 1]
    bigind_file_start = file_ibeam_start
    bigind_file_end = file_ibeam_end

    @save "file_beam_inds.jld2" bigind_file_start nbeams bigind_file_end
    return file_ibeam_start, file_ibeam_end
end

if !isfile("lidar_dt.jld2")
    error("Expected lidar_dt.jld2 to exist. Run save_lidar_dt.jl first.")
end

lidardata = load("lidar_dt.jld2")
dtime = lidardata["dtime"]
ist = Int.(lidardata["ist"])
ien = Int.(lidardata["ien"])

lidarstemdir = "./data"
ncdir = joinpath(lidarstemdir, "netcdf_stare")
ncfiles = sort(filter(f -> startswith(f, "Stare_116_") && endswith(f, ".nc"), readdir(ncdir)))
files = joinpath.(ncdir, ncfiles)

file_ibeam_start, file_ibeam_end = load_file_beam_inds(files)

idx = lidar_index.build_lidar_index(dtime, ist, ien, files, file_ibeam_start, file_ibeam_end)
@save "lidar_index.jld2" idx

println("Saved lidar_index.jld2 with ", length(idx.ibeam_start), " chunks across ", length(idx.filenames), " files.")

# @load "lidar_index.jld2" idx