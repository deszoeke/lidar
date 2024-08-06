# MWE for finding and sorting start and end of stares

using Dates
include("read_lidar.jl")
using .read_lidar

fullfiles = joinpath.(pwd(),"data","20240531",files)
ta, hdr, nbeams = read_lidar.read_streamlinexr_beam_timeangles(fullfiles)
ien = findall( diff(ta[:time]) .> 0.01 ) # gap of more than a few tens of seconds
ist = ien .+ 1 # start of the next
# ii increasing-order vector of all start and end indices
ii = sort( permutedims([ien ist])[:] ) # this should sort them, but sort() for good measure
isend = in.(ii, Ref(ien)) # Boolean true for interval end points of ii