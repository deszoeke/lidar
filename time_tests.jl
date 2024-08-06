# find and sorting start and end of stare chunks

using Revise
using Dates
include("read_lidar.jl")
using .read_lidar


# Request chunk subscripts and start and end times in any time window.

# get pickets marking start and end of chunks for one day
datapath = joinpath.(pwd(),"data","20240531")
files = filter(startswith("Stare"), readdir(datapath))
fullfiles = joinpath.(datapath, files)
ta, hdr, nbeams = read_lidar.read_streamlinexr_beam_timeangles(fullfiles)

# identify chunks from gaps
ien = findall( diff(ta[:time]) .> 0.01 ) # gap of more than a few tens of seconds
ist = ien .+ 1 # start of the next chunk

# ii increasing-order vector of all start and end indices
ii = sort( permutedims([ien ist])[:] ) # ordered; this should sort them, but sort() for good measure
jj = LinearIndices(ii) # jj indexes ii
# ii alternates end start-end ... start ii[jjen] == ien; ii[jjst] == ist
jjen = jj[1:2:end]
jjst = jj[2:2:end] # even starts of next chunks follow previous odd ends
# NOTE: a view can change the parity of starts and ends

# functions to check if an index i is an end point; or a start point
isend(i) = in(i, Ref(ien)) # Boolean true for interval chunk end points of ii
isstart(i) = in(i, Ref(ist)) # Boolean true for interval chunk start points of ii

# end index corresponding to a start indes
# get the index j such that i = ii[j]
thisj(i, ii=ii) = findfirst(ii.==i)
# get the end index ien corresponding to the start index ist of the same chunk
chunken(ist, ii=ii) = ii[thisj(ist, ii) + 1]
# get the start index ist corresp to the end index ien of the same chunk
chunkst(ien, ii=ii) = ii[thisj(ien, ii) - 1]

# index the next chunk
nextchunki(i, ii=ii) = ii[thisj(i, ii) + 2]
# index the previous chunk
prevchunki(i, ii=ii) = ii[thisj(i, ii) - 2]

"""
Find chunks within query start and end times.
Inclusive: include chunks for which query falls in the middle.
"""
function query_start_end_chunks(qst, qen; ien=ien, ist=ist)
    ien1 = ien[findfirst(x-> x>=qst , ien)] # end index of first queried chunk     # do in this order
    ist1 = ist[findlast( x-> x <ien1, ist)] # start index of first queried chunk
    lastist = ist[findlast( x-> x<=qen   , ist)] # start index of last queried chunk
    lastien = ien[findfirst(x-> x>lastist, ien)] # end index of last queried chunk
    ist1,ien1, lastist,lastien
end

(ist1,ien1, lastist,lastien) = query_start_end_chunks(qst, qen; ien=ien, ist=ist)
# pickets of queried subset
iisub = ii[ ist1 .>= ii .>= lastien ] # start-end ... start-end ; no orphans
# iisub parity switched from ii
istsub = iisub[1:2:end]
iensub = iisub[2:2:end] # even ends follow odd starts of same chunk

##

# Which file(s) is a chunk in?
# time index hierarchy:
#
# beam
#   chunk [start, end]
#       file [start, end]

# Collect data start indices times of all files in ibeam1
# number of beams in each file nbeam
# so last beam index in a file is ibeam1+nbeam-1
get_start_fileidx(tidx) = findlast(t-> ibeam1<=t, tidx)
get_end_fileidx(tidx) = findfirst(t-> ibeam1+nbeam-1>=t, tidx)
get_stare_files(tidx, files) = files[get_fileidx(tidx)]

