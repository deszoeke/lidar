# Construct and save a time vector of datetimes for all vertically staring lidar beams
# Reads from data/netcdf_stare/*.nc (not .hpl) so April 30 coverage is complete.

# preamble
using Revise
using Pkg; Pkg.activate(".")

using Dates
using JLD2
using NCDatasets

include("./read_lidar.jl")
using .read_lidar
using .read_lidar.stare

lidarstemdir = "./data"
ncdir = joinpath(lidarstemdir, "netcdf_stare")

# Read DateTime time vectors from all NC stare files in sorted order.
# NC filenames have the correct calendar date, and each file's :time is already
# a proper DateTime — no decimal-hour arithmetic or day-wrap heuristics needed.
ncfiles = sort(filter(f -> startswith(f, "Stare_116_") && endswith(f, ".nc"), readdir(ncdir)))
ff = joinpath.(ncdir, ncfiles)
println("Reading $(length(ff)) NetCDF stare files...")

dtime = DateTime[]
nbeams_per_file = Int[]
for (i, f) in enumerate(ff)
    beams = read_lidar.read_stare_hourly_netcdf(f; vars=[:time])
    t = beams[:time]
    append!(dtime, t)
    push!(nbeams_per_file, length(t))
    i % 24 == 0 && println("  $(i)/$(length(ff)): $(basename(f)) → $(length(t)) beams")
end
dtime = convert(Vector{DateTime}, dtime)
println("Total beams: $(length(dtime))  from $(dtime[1]) to $(dtime[end])")

# indices of gaps to find stare chunks
# ist/ien are FULL chunk arrays (length = nchunks):
#   ist[ic] = first global beam of chunk ic,  ien[ic] = last global beam of chunk ic
function stare.all_gaps(dt::Vector{DateTime})
    gap_ends   = findall( diff(dt) .> Second(36) )   # last beam before each gap
    gap_starts = gap_ends .+ 1                        # first beam after each gap
    ist = [1;           gap_starts]                   # chunk starts: beam 1, then after each gap
    ien = [gap_ends;    length(dt)]                   # chunk ends:   before each gap, then last beam
    return ist, ien
end
ist, ien = all_gaps(dtime)

@save "lidar_dt.jld2" dtime ist ien
println("Saved lidar_dt.jld2: $(length(dtime)) beams, $(length(ist)) chunks")