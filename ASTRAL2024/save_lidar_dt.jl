# Construct and save a time vector of datetimes for all vertically staring lidar beams

# preamble
using Revise
using Pkg; Pkg.activate(".")

using Dates
using JLD2 

include("./read_lidar.jl")
using .read_lidar
using .read_lidar.stare
# using .read_vecnav: read_vecnav_dict
# import .chunks
# include("./timing_lidar.jl")
# using .timing_lidar

using PyPlot
using PyCall
using PyCall: PyObject

# PyObject method interprets Array{Union{T,Missing}} as a
# numpy masked array.
# This allows for plotting with missing values.
function PyObject(a::Array{Union{T,Missing},N}) where {T,N}
    numpy_ma = PyCall.pyimport("numpy").ma
    pycall(numpy_ma.array, Any, coalesce.(a,zero(T)), mask=ismissing.(a))
end

# get all times to make a main list of stare chunks, regardless of file boundaries
starefiles = filter(startswith("Stare_116_"), readdir(joinpath(lidarstemdir, "all")))
ff = joinpath.(lidarstemdir, "all", starefiles)
ta, _,_ = read_lidar.read_streamlinexr_beam_timeangles(ff) # slow

day0 = floor(ta[:start_time][1], Day)
dtime = day0 .+ Millisecond.(round.(Int64, 3600_000 .* ta[:time]))
for j in findall(diff(ta[:time]) .< -1.0)
    dtime[(j+1):end] .+= Day(1) # increment days
end

# indices of gaps to find stare chunks
function stare.all_gaps(dt::Vector{DateTime})
    ien = findall( diff(dt) .> Second(36) )
    ist = ien .+ 1
    return ien, ist
end
ien, ist = all_gaps(dtime)

@save "lidar_dt.jld2" dtime ist ien

load