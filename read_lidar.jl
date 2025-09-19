# using Revise
# using Pkg; Pkg.activate("/Users/deszoeks/Projects/ASTRAL/lidar")

# contains code for several modules:
# read_lidar
# read_lidar.stares
# read_vecnav
# chunks

module read_lidar

using Dates
using NCDatasets
using Interpolations
using JLD2
# using MAT
# using PyPlot

export read_streamlinexr_stare #, read_streamlinexr_head
export get_daily_meanuv, read_daily_Vn
export read_streamlinexr_beam_timeangles

export read_lidar_chunks # used
export get_all_file_start_end_idxs # exported to test

### utilites
pd = permutedims
m2n(x) = ismissing(x) ? NaN : x

### read stares functions
# functions for reading hpl files
# # lidardir = "./data/lidar"
lidardir = "./data"

"count the headlines before a line that starts with ****"
function numheaderlines(file)
    n = 0 
    open(file) do f
        while !startswith( readline(f), "****" )
            n += 1
        end
    end
    return n + 1
end

"""
header = read_streamlinexr_head(file_path)
Read the Halo Photonics StreamLineXR .hpl file header.
"""
function read_streamlinexr_head(file_path)
    # put the header in a Dict
    h = Dict(
        :nlines => countlines(file_path),
        :nheaderlines => 
        # variables defined outside the do block, defaults will be overwritten.
        :ngates => Int32(104),
        :gatelength => 30.0,
        :nrays => Int32(1), # or 40
        :start_time => DateTime(2024,5,2,10,0,0),
        :vel_resolution => 0.0380,
        :headr1 => split("time Azimuth Elevation Pitch Roll"),
        :units1 => split("hours degrees degrees degrees degrees"),
        :headr2 => ["Range Gate", "Doppler", "Intensity", "Beta"],
        :units2 => ["none", "(m/s)", "(SNR + 1)", "(m-1 sr-1)"]
    )

    # open and read the file header, updating the header dict a
    open(file_path) do file
        for _ in 1:2 # skip lines
            readline(file)
        end

        # read number of gates from line 3
        h[:ngates] = parse(Int32, last(split(readline(file)))) # 104
        # read gate length from line 4
        h[:gatelength] = parse(Float64, last(split(readline(file))))
        for _ in 5:6 # skip
            readline(file)
        end

        # read number of rays from line 7
        h[:nrays] = parse(Int32, last(split(readline(file)))) # 1 for Stare or 40 for User
        for _ in 8:9 # skip
            readline(file)
        end

        # read start time from line 10
        tl = split(readline(file))
        h[:start_time] = DateTime(tl[end-1]*tl[end], "yyyymmddHH:MM:SS.ss")
        # read velocity resolution (m/s) line 11
        h[:vel_resolution] = parse(Float64, last(split(readline(file))))
        readline(file) # skip line 12
        # read header 1 from line 13
        hdrln = split(readline(file))[5:end]
        h[:headr1] = hdrln[1:2:end]
        h[:units1] = hdrln[2:2:end]
        readline(file) # skip line 14
        # read header 2 from line 15
        hdrln = split(readline(file))[4:end] # but don't parse it
        h[:headr2] = ["Range Gate", "Doppler", "Intensity", "Beta"]
        h[:units2] = ["none", "(m/s)", "(SNR + 1)", "(m-1 sr-1)"]
        # for _ in 16:17 # skip
        #     readline(file)
        # end
    end

    return h
end

"""
modifying read_streamlinexr_stare!(file_path, header, beams)
Read data and fill in the beams for a single file.
"""
function read_streamlinexr_stare!(file_path, h, beams, nheaderlines=17; nbeams0=0, startat=1, endat=0)
    # nbeams0 is the number of beams alread read into the Dict beams
    #   writing to the Dict fills in from there.
    # nbeamsmax = fileen[fj]


    # use header information in h
    nlines = h[:nlines]
    ngates = h[:ngates]

    # beams could be rays or times
    nbeamsmax = round(Int, (nlines-nheaderlines) / (1+ngates)) # = nrays*ntimes # total number available
    # but nbeams may be reduced by startat, endat
    endat = mod(endat-1, nbeamsmax) + 1 # clobbers name but does not write to original argument
    nbeams = min(endat - startat + 1, nbeamsmax) # actual number of beams requested, or total number available

    # allocates for each file; this is not too much to affect perfomance
    beam_timeangles = zeros(Float64, (nbeams, 5))
    beam_velrad = zeros(Float64, nbeams, ngates, 4)

    # for User wind profiles beam <--> VAD ray
    # for Stare beam <--> time

    # open and read the file
    open(file_path) do file
        for _ in 1:nheaderlines # skip header lines
            readline(file)
        end
        for _ in 1:( (1+ngates) * (startat-1) ) # skip beams before startat
            readline(file)
        end

        # now read data # nbeams is already limited by endat-startat+1
        for ibeam = 1:nbeams
            # beam described by a batch of 1+ngates lines
            # Read the beam parameter line
            line = readline(file)
            try
                beam_timeangles[ibeam,:] .= parse.(Float64, split(line))
            catch
                @show line
            end
            # Read each gate in the beam
            for igate = 1:ngates
                line = readline(file)
                beam_velrad[ibeam, igate,:] .= parse.(Float64, split(line))
            end
        end
    end # close the file

    # parse the variables into the dict beams
    # by beam
    # offset bb by beam0, indices already read into Dict beam
    bb = nbeams0 .+ (1:nbeams) # for the output array dimension nbeams
    beams[:time     ][bb] .= beam_timeangles[:,1] # decimal hours
    beams[:azimuth  ][bb] .= beam_timeangles[:,2] # degrees
    beams[:elevangle][bb] .= beam_timeangles[:,3] # degrees
    beams[:pitch    ][bb] .= beam_timeangles[:,4]
    beams[:roll     ][bb] .= beam_timeangles[:,5]
    # by gate
    beams[:height   ] .= (beam_velrad[1,:,1].+0.5) .* h[:gatelength] # center of gate, assumes same for all beams

    # dependent variables (beam, gate)
    beams[:dopplervel][bb,:] .= beam_velrad[:,:,2] # m/s
    beams[:intensity ][bb,:] .= beam_velrad[:,:,3] # SNR + 1
    beams[:beta      ][bb,:] .= beam_velrad[:,:,4] # m-1 sr-1  backscatter?
end

"""
beams, h = read_streamlinexr_stare(file_path)
Read the hourly Photonics StreamLineXR Stare_... .hpl file data and header.
"""
# function read_streamlinexr_stare( d::Date )
#     # index by directory contents
#     files = readdir(joinpath(lidardir,datestamp,starefile)
#     read_streamlinexr_stare( file_path )
# end
# !! Must specify hourly files.
function read_streamlinexr_stare( dt::DateTime )
    yyyymmdd    = Dates.format( dt, dateformat"yyyymmdd" )
    yyyymmdd_HH = Dates.format( dt, dateformat"yyyymmdd_HH" )
    starefile = "Stare_116_$(yyyymmdd_HH).hpl"
    file_path = joinpath(lidardir,yyyymmdd,starefile)
    read_streamlinexr_stare( file_path )
end

"read a single Stare file"
function read_streamlinexr_stare(file_path::AbstractString, nheaderlines=17; startat_=1, endat_=0)
    # use header information in h 
    h = read_streamlinexr_head(file_path)
    nlines = h[:nlines]
    ngates = h[:ngates]

    # for single-file
    # beams contain data from rays at particular times
    nbeamsmax = round(Int, (nlines-nheaderlines) / (1+ngates)) # = nrays*ntimes # total number available
    # nbeams may be reduced by startat, endat
    startat = max(1, min( mod(startat_-1, nbeamsmax) + 1, nbeamsmax )) # mod provides max min limits anyway
    endat = mod(endat_-1, nbeamsmax) + 1 # 0 goes to nbeamsmax (end)
    nbeams = min(endat - startat + 1, nbeamsmax) # actual number of beams requested, or total number available

    # initialize a beams Dict with exactly nbeams
    beams = Dict(
        :time      => Vector{Union{Float32,Missing}}(missing, nbeams), # decimal hours
        :azimuth   => Vector{Union{Float32,Missing}}(missing, nbeams), # degrees
        :elevangle => Vector{Union{Float32,Missing}}(missing, nbeams),
        :pitch     => Vector{Union{Float32,Missing}}(missing, nbeams),
        :roll      => Vector{Union{Float32,Missing}}(missing, nbeams),

        :height    => Vector{Union{Float32,Missing}}(missing, ngates), # center of gate

        # dependent variables (beam, gate)
        :dopplervel => Matrix{Union{Float32,Missing}}(missing, nbeams,ngates), # m/s
        :intensity  => Matrix{Union{Float32,Missing}}(missing, nbeams,ngates), # SNR + 1
        :beta       => Matrix{Union{Float32,Missing}}(missing, nbeams,ngates) # m-1 sr-1  backscatter?
        )

    # read file and fill beams with data
    # read_streamlinexr_stare!(file_path, h, beams, nheaderlines)
    read_streamlinexr_stare!(file_path, h, beams, nheaderlines; nbeams0=0, startat=startat, endat=endat)
    return beams, h
end

"read multiple Stare files"
function read_streamlinexr_stare(file_path::AbstractVector{T}, fi1, filast, bigstartat=fi1[1], bigendat=last(filast), nheaderlines=17) where {T<:AbstractString}
    # loop over a vector of files.
    # bigstartat (default fi1[1]) is the first requested start index of all the data in the files
    # bigendat (default last(filast)) is the last requested end index of all the data in the files.

    # set up
    mask = ( filast.>=bigstartat .&& fi1.<=bigendat )
    file_path = file_path[ mask ] # trim just keeping the files we need
    fi1       = fi1[       mask ]
    filast    = filast[    mask ]

    nfiles = length(file_path)
    nbeams  = zeros(Int32, nfiles)  # number of beams in each file
    h = read_streamlinexr_head(file_path[1]) # use header information in h
    ngates = h[:ngates]

    # number of beams total already known outside and can be calculated
    nb = bigendat - bigstartat + 1 # total number of beams in requested stares
    # initialize a beams Dict
    beams = Dict(
        :time      => Vector{Union{Float32,Missing}}(missing, nb), # decimal hours
        :azimuth   => Vector{Union{Float32,Missing}}(missing, nb), # degrees
        :elevangle => Vector{Union{Float32,Missing}}(missing, nb),
        :pitch     => Vector{Union{Float32,Missing}}(missing, nb),
        :roll      => Vector{Union{Float32,Missing}}(missing, nb),

        :height    => Vector{Union{Float32,Missing}}(missing, ngates), # center of gate

        # dependent variables (beam, gate)
        :dopplervel => Matrix{Union{Float32,Missing}}(missing, nb,ngates), # m/s
        :intensity  => Matrix{Union{Float32,Missing}}(missing, nb,ngates), # SNR + 1
        :beta       => Matrix{Union{Float32,Missing}}(missing, nb,ngates) # m-1 sr-1  backscatter?
        )

    # read each file and write data in Dict beams
    nbeams0 = 0
    for (fi, file) in enumerate(file_path)
        h = read_streamlinexr_head(file) # reread header for each file
        nlines = h[:nlines]
        nbeamsmax = round(Int, (nlines - nheaderlines) / (1+ngates)) # number of beams for each file
        # set startat, endat if in first and/or last file, otherwise 1, nbeamsmax
        if fi == 1
            startat = bigstartat - fi1[fi] + 1 # index referenced to file
        else
            startat = 1
        end
        if fi == nfiles
            endat = bigendat - fi1[fi] + 1
        else
            endat = nbeamsmax
        end
        # read and update beams[keys][nbeams0 .+ (1:nbeams[i])]
        read_streamlinexr_stare!(file, h, beams, nheaderlines; nbeams0=nbeams0, startat=startat, endat=endat) 
        nbeams0 += nbeams[fi] # increment number of beams now read
    end
    return beams, h, nbeams0
end

"""
Read multiple whole files.
provide no refined query indices.
"""
function read_streamlinexr_stare(file_path::AbstractVector{T}, nheaderlines=17) where {T<:AbstractString}
    # loop over a vector of files.
    
    # use header information in h, count lines in all files
    nfiles = length(file_path)
    nbeams = zeros(Int32, nfiles)
    ngates = -1
    # read number of lines for each file
    for (i,file) in enumerate(file_path)
        h = read_streamlinexr_head(file)
        nlines = h[:nlines]
        ngates = h[:ngates]
        # beams could be rays or times
        nbeams[i] = round(Int, (nlines - nheaderlines) / (1+ngates)) # number of beams for each file
    end

    # initialize a beams Dict
    nb = sum(nbeams) # total number of beams in all files
    beams = Dict(
        :time      => Vector{Union{Float32,Missing}}(missing, nb), # decimal hours
        :azimuth   => Vector{Union{Float32,Missing}}(missing, nb), # degrees
        :elevangle => Vector{Union{Float32,Missing}}(missing, nb),
        :pitch     => Vector{Union{Float32,Missing}}(missing, nb),
        :roll      => Vector{Union{Float32,Missing}}(missing, nb),

        :height    => Vector{Union{Float32,Missing}}(missing, ngates), # center of gate

        # dependent variables (beam, gate)
        :dopplervel => Matrix{Union{Float32,Missing}}(missing, nb,ngates), # m/s
        :intensity  => Matrix{Union{Float32,Missing}}(missing, nb,ngates), # SNR + 1
        :beta       => Matrix{Union{Float32,Missing}}(missing, nb,ngates) # m-1 sr-1  backscatter?
        )

    # read each file and fill beams with data
    nbeams0 = 0
    for (i,file) in enumerate(file_path)
        h = read_streamlinexr_head(file) # reread header for each file
        read_streamlinexr_stare!(file, h, beams, nheaderlines; nbeams0=nbeams0) # updates beams[keys][nbeams0 .+ (1:nbeams[i])]
        nbeams0 += nbeams[i] # number of beams now read
    end
    return beams, h, nbeams0
end

"update beam time and angle data read from a single file"
function read_streamlinexr_beam_timeangles!(file_path, h, beams, nheaderlines=17; nbeams0=0)
    # use header information in h
    nlines = h[:nlines]
    ngates = h[:ngates]

    # beams could be rays or times
    nbeams = round(Int, (nlines-nheaderlines) / (1+ngates)) # = nrays*ntimes
    beam_timeangles = zeros(nbeams, 5)

    # for User wind profiles beam <--> VAD ray
    # for Stare beam <--> time

    # open and read the file
    open(file_path) do file
        for _ in 1:nheaderlines # skip header lines
            readline(file)
        end

        # now read parameter data for each beam
        for ibeam = 1:nbeams
            line = readline(file)
            beam_timeangles[ibeam,:] .= parse.(Float64, split(line))
            # skip reading the beam velocity and backscatter data
            [readline(file) for i in 1:ngates]
        end
    end # close the file

    # put the beam parameter variables into the dict beams
    bb = nbeams0 .+ (1:nbeams) # for the output array dimension nbeams
    beams[:time     ][bb] .= beam_timeangles[:,1] # decimal hours
    beams[:azimuth  ][bb] .= beam_timeangles[:,2] # degrees
    beams[:elevangle][bb] .= beam_timeangles[:,3] # degrees
    beams[:pitch    ][bb] .= beam_timeangles[:,4]
    beams[:roll     ][bb] .= beam_timeangles[:,5]
end


"read multiple file time and angle data"
function read_streamlinexr_beam_timeangles(file_path::AbstractVector{T}, nheaderlines=17) where {T<:AbstractString}
    # loop over a vector of files.
    
    # use header information in h, count lines in all files
    nfiles = length(file_path)
    nbeams = zeros(Int32, nfiles)
    ngates = -1
    # read number of lines for each file
    for (i,file) in enumerate(file_path)
        h = read_streamlinexr_head(file)
        nlines = h[:nlines]
        ngates = h[:ngates]
        # beams could be rays or times
        nbeams[i] = round(Int, (nlines - nheaderlines) / (1+ngates)) # number of beams for each file
    end

    # initialize a beams Dict
    nb = sum(nbeams) # total number of beams in all files
    beams = Dict(
        :time      => Vector{Union{Float32,Missing}}(missing, nb), # decimal hours
        :azimuth   => Vector{Union{Float32,Missing}}(missing, nb), # degrees
        :elevangle => Vector{Union{Float32,Missing}}(missing, nb),
        :pitch     => Vector{Union{Float32,Missing}}(missing, nb),
        :roll      => Vector{Union{Float32,Missing}}(missing, nb),
        :ibeam1    =>  Vector{Union{Int32,   Missing}}(missing, nfiles), # index of first beam in file
        :nbeams     => Vector{Union{Int32,   Missing}}(missing, nfiles),  # number of beams in each file
        :start_time => Vector{Union{DateTime,Missing}}(missing, nfiles) # include base time
        )

    # read each file and fill beams with data
    nbeams0 = 0
    h = read_streamlinexr_head(file_path[1]) # first header
    for (fi, file) in enumerate(file_path)
        h = read_streamlinexr_head(file) # reread header for each file
        read_streamlinexr_beam_timeangles!(file, h, beams, nheaderlines; nbeams0=nbeams0) # updates beams[keys][nbeams0 .+ (1:nbeams[i])]
        beams[:ibeam1][fi] = nbeams0 + 1 # index of first beam in file
        beams[:nbeams][fi] = nbeams[fi]   # number of beams in this file
        beams[:start_time][fi] = h[:start_time]
        nbeams0 += nbeams[fi] # number of beams now read
    end
    return beams, h, nbeams0
end

# just reading time data is not actually quicker than timeangles
"update beam time data read from a single file"
function read_streamlinexr_beam_times!(file_path, h, times, nheaderlines=17; nbeams0=0)
    # use header information in h
    nlines = h[:nlines]
    ngates = h[:ngates]

    # beams could be rays or times
    nbeams = round(Int, (nlines-nheaderlines) / (1+ngates)) # = nrays*ntimes
    # beam_timeangles = zeros(nbeams, 5)

    # open and read the file
    open(file_path) do file
        for _ in 1:nheaderlines # skip header lines
            readline(file)
        end

        # now read parameter data for each beam
        for ibeam = 1:nbeams
            line = readline(file)
            times[nbeams0+ibeam, :] .= parse.(Float64, split(line)[1]) # just read time
            # skip reading the beam velocity and backscatter data
            [readline(file) for i in 1:ngates]
        end
    end # close the file

    # put the beam parameter variables into the dict beams
    # bb = nbeams0 .+ (1:nbeams) # for the output array dimension nbeams
    # beams[:time     ][bb] .= beam_timeangles[:,1] # decimal hours
    # beams[:azimuth  ][bb] .= beam_timeangles[:,2] # degrees
    # beams[:elevangle][bb] .= beam_timeangles[:,3] # degrees
    # beams[:pitch    ][bb] .= beam_timeangles[:,4]
    # beams[:roll     ][bb] .= beam_timeangles[:,5]
end


"read multiple file time data"
function read_streamlinexr_beam_times(file_path::AbstractVector{T}, nheaderlines=17) where {T<:AbstractString}
    # loop over a vector of files.
    
    # use header information in h, count lines in all files
    nfiles = length(file_path)
    nbeams = zeros(Int32, nfiles)
    ngates = -1
    # read number of lines for each file
    for (i,file) in enumerate(file_path)
        h = read_streamlinexr_head(file)
        nlines = h[:nlines]
        ngates = h[:ngates]
        # beams could be rays or times
        nbeams[i] = round(Int, (nlines - nheaderlines) / (1+ngates)) # number of beams for each file
    end

    # initialize a beams Dict
    nb = sum(nbeams) # total number of beams in all files
    # beams = Dict(
    #     :time      => Vector{Union{Float32,Missing}}(missing, nb), # decimal hours
    #     :azimuth   => Vector{Union{Float32,Missing}}(missing, nb), # degrees
    #     :elevangle => Vector{Union{Float32,Missing}}(missing, nb),
    #     :pitch     => Vector{Union{Float32,Missing}}(missing, nb),
    #     :roll      => Vector{Union{Float32,Missing}}(missing, nb),
    #     :ibeam1    => Vector{Union{Int32,Missing}}(missing, nfiles), # index of first beam in file
    #     :nbeams    => Vector{Union{Int32,Missing}}(missing, nfiles)  # number of beams in each file
    #     )
    times  = Vector{Union{Float32,Missing}}(missing, nb)    # decimal hours
    ibeam1 = Vector{Union{Int32,Missing}}(missing, nfiles)  # index of first beam in file
    # read each file and fill beams with time and index data
    nbeams0 = 0
    h = read_streamlinexr_head(file_path[1]) # first header
    for (fi, file) in enumerate( file_path )
        h = read_streamlinexr_head(file) # reread header for each file
        read_streamlinexr_beam_times!(file, h, times, nheaderlines; nbeams0=nbeams0) # updates beams[keys][nbeams0 .+ (1:nbeams[i])]
        ibeam1[fi] = nbeams0 + 1 # index of first beam in file
        nbeams0 += nbeams[fi]    # number of total beams now read
    end
    return times, ibeam1, nbeams, h  # note changed outputs
end

### read mean wind functions
# uvdir = "./data/lidar/netcdf/"
uvdir = "./data/netcdf/"

"get daily mean-wind netcdf dataset"
function get_daily_meanuv( dt::Union{Date, DateTime} )
    get_daily_meanuv( Dates.format(Date(dt), dateformat"yyyymmdd") )
end
function get_daily_meanuv( yyyymmdd::AbstractString )
    Dataset(joinpath(uvdir,"ekamsat_lidar_uv_$(yyyymmdd).nc"))
end

### read VectorNav functions
# JLD2 files created by vectornav.ipynb

# Vndir = "./data/table/" # uses symbolic link ./data in cwd
# no daily leg2 jld2 files; use 

"read daily JLD2 as a Dict"
function read_daily_Vn( dt::Union{Date, DateTime} )
    Vndir = dt < DateTime(2024,5,31) ? "./data/table/leg1" : "./data/table/leg2"
    yyyymmdd = Dates.format(Date(dt), dateformat"yyyymmdd")
    read_daily_Vn( joinpath(Vndir, "VectorNavTable_$(yyyymmdd).jld2") )
end
function read_daily_Vn( filename::AbstractString )
    Vn = JLD2.load( filename )
end

## utility functions to determine which file a chunk in and its indices in that file

# !MOVE utilities for requesting, reading, and subsetting data from files
# TO read_lidar.jl
# include("read_lidar.jl")

# internal functions:
# get_stare_files, thisfile_idx, idx_beams_in_file, thisfile_start_end, thisfile_name
# export read_lidar_chunks # used
# export get_all_file_start_end_idxs # exported to test

# time index hierarchy:
#     beam
#   stare chunk [start, end]
# file [start, end]

# Collect data start indices times of all files in ibeam1
# number of beams in each file nbeam
# so last beam index in a file is ibeam1+nbeam-1
# get_start_fileidx(tidx) = findlast(t-> ibeam1<=t, tidx)
# get_end_fileidx(tidx) = findfirst(t-> ibeam1+nbeam-1>=t, tidx)
get_stare_files(tidx, files) = files[get_fileidx(tidx)]

# vectors of all file start and end indices
function get_all_file_start_end_idxs(ta::Dict)
    fi1 = ta[:ibeam1][:]
    filast = @. fi1 + ta[:nbeams] - 1
    return fi1, filast
end
function get_all_file_start_end_idxs(file_paths) 
    ta = read_streamlinexr_beam_timeangles(file_path)
    get_all_file_start_end_idxs(ta)
end

"index of the file for which the index x falls between fi1 and filast"
thisfile_idx(x, fi1, filast) = findfirst( fi1.<= x .<=filast)

# this method superseded by one that uses fi1, filast in calling scope
idx_beams_in_file(i, ibeam1, fi1, filast) = i - ibeam1[thisfile_idx(i, fi1, filast)] + 1

function thisfile_start_end(x, fi1, filast)
    j = thisfile_idx(x, fi1, filast)
    return fi1[j], filast[j]
end

"get name of a file corresponding to a time or index"
thisfile_name(x::Integer, fi1, filast, files=fullfiles) = files[thisfile_idx(x, fi1, filast)]
function thisfile_name(dtx::DateTime, dt::Vector{DateTime}, files=fullfiles) 
    fi1, filast = get_all_file_start_end_idxs(files)
    thisfile_name( findlast(dt .<= dtx), fi1, filast, files)
end

"read chunks inclusive between [firstdtq lastdq] from vector of files"
function read_lidar_chunks(files, firstdtq, lastdtq)
    # read all times
    ta, hdr, nbeams = read_lidar.read_streamlinexr_beam_timeangles(files) # 5 s per day

    # TODO append day hour to the decimal hour...
    hdr[:start_time]
    # finding indices
    ien, ist = all_gaps(ta[:time]) # identify indices of chunks from gaps
    ii, jj = all_start_end_indices(ien, ist)

    # Request some beam times. Low-level query takes indices.
    # # find indices qst, qen for the request
    # dt = @. DateTime(2024,5,31) + Millisecond(round(Int64, ta[:time] * 3_600_000)) # ta[:time] is decimal hour
    # # problem: date not defined in ta
    # qst = findfirst(dt .>= DateTime(2024,5,31,5))
    # qen = findlast( dt .<  DateTime(2024,5,31,10))
    qst = findfirst(dt .>= firstdtq)
    qen = findlast( dt .<  lastdtq)

    # request the first and last chunks
    (ist1,ien1, lastist,lastien) = query_start_end_chunks(qst, qen; ien=ien, ist=ist)
    # request start and end of all the chunks between ist1 and lastien
    istsub, iensub = all_chunks(ist1, lastien)

    # sort the chunks by file

    # Subset the data within the file to the requested chunks
    # could be done by skipping to start_beam, and stopping at stop_beam.
    ibeam1 = ta[:ibeam1]
    # not needed:
    ## idx_beams_in_file(i, ibeam1) = i - ibeam1[thisfile_idx(i, ist1, lastien)] + 1
    # start_beam = idx_beams_in_file(ist1, ibeam1, ist1, lastien) # default = 1
    # stop_beam = idx_beams_in_file(lastien, ibeam1, ist1, lastien)      # default = ta[:nbeams]
    
    # only ever truncate the data read by a file because of a
    # intentional scientific user-level query

    # @show ist1, lastien
    # @show [fi1 filast]

    # read the beam data between ist1 and lastien
    beams, h, nbeams0 = read_lidar.read_streamlinexr_stare(files, fi1, filast, ist1, lastien)
    return beams, h, ien, ist
end

# stare functions to find and use start and end indices of stare chunks

module stare

using Dates

# used in read_lidar_chunks and tests
export hour_beams
export all_gaps
export all_start_end_indices
export query_start_end_chunks
export all_chunks

# Request chunk subscripts and start and end times in any time window.
"return time as a decimal hour from the start time"
function hour_beams(time, ibeam1, start_time)
    # time is decimal hour of the day
    # ibeam1 is the index of time corresponding to start_time
    # start_time is the DateTime when a new file starts
    
    # expand time axis to include day and decimal hour
    ibeam_expand = [ findlast( ibeam1 .<= i ) for i in eachindex(time) ]
    dtbeams = @. floor(start_time, Day)[ibeam_expand] 
    hourbeams = @. 24*Dates.value(Day(dtbeams-dtbeams[1])) + time
end
"hour_beams(timeangles::Dict)"
hour_beams(ta::Dict) = hour_beams(ta[:time], ta[:ibeam1], ta[:start_time])

# identify indices of chunks from gaps
function all_gaps(t)
    ien = findall( diff(t) .> 0.01 ) # gap of more than a few tens of seconds
    ist = ien .+ 1 # start of the next chunk
    return ien, ist
end

function all_gaps(dt::Vector{DateTime})
    ien = findall( diff(dt) .> Second(36) )
    ist = ien .+ 1
    return ien, ist
end

function all_start_end_indices(ien, ist)
    # ii increasing-order vector of all start and end indices
    ii = sort( permutedims([ien ist])[:] ) # ordered; this should sort them, but sort() for good measure
    jj = LinearIndices(ii) # jj indexes ii
    # ii alternates end start-end ... start ii[jjen] == ien; ii[jjst] == ist
    jjen = jj[1:2:end]
    jjst = jj[2:2:end] # even starts of next chunks follow previous odd ends
    # NOTE: a view can change the parity of starts and ends
    return ii,jj, jjen,jjst
end
# ii is the big array of start and end indices

# functions to check if an index i is an end point; or a start point
isend(i) = in(i, Ref(ien)) # Boolean true for interval chunk end points of ii
isstart(i) = in(i, Ref(ist)) # Boolean true for interval chunk start points of ii

# NOTE ii default is computed by all_start_end_indices, but not available in the lexical scope
# so ii MUST be provided in the call

# end index corresponding to a start indes
# get the index j such that i = ii[j]
thisj(i, ii) = findfirst(ii.==i)
# get the end index ien corresponding to the start index ist of the same chunk
chunken(ist, ii) = ii[thisj(ist, ii) + 1]
# get the start index ist corresp to the end index ien of the same chunk
chunkst(ien, ii) = ii[thisj(ien, ii) - 1]

# index the next chunk
nextchunki(i, ii) = ii[thisj(i, ii) + 2]
# index the previous chunk
prevchunki(i, ii) = ii[thisj(i, ii) - 2]

# Find chunks within query start and end times.
"""
Request start and end chunks based on time indices.
Inclusive: include chunks for which query falls in the middle.
"""
function query_start_end_chunks(qst, qen; ien=ien, ist=ist)
    ien1 = ien[findfirst(x-> x>=qst , ien)] # end index of first queried chunk     # do in this order
    ist1 = ist[findlast( x-> x <ien1, ist)] # start index of first queried chunk
    lastist = ist[findlast( x-> x<=qen   , ist)] # start index of last queried chunk
    lastien = ien[findfirst(x-> x>lastist, ien)] # end index of last queried chunk
    ist1,ien1, lastist,lastien
end

"""
request all start indices and all end indices between 
index of start of first chunk and end of last chunk lastien
"""
function all_chunks(ist1, lastien, ii)
    # pickets of queried subset
    iisub = ii[ findall(x-> ist1<=x<=lastien, ii) ] # start-end ... start-end ; no orphans

    istsub = iisub[1:2:end] # iisub parity switched from ii
    iensub = iisub[2:2:end] # even ends follow odd starts of same chunk
    return istsub, iensub # start-end pairs
end

end # module stare

end # module read_lidar


module read_vecnav
# functions to read the VectorNav file # MOVE TO read_lidar

using Dates
using JLD2

#using .timing_lidar # mightcould wanna use this

export read_vecnav_dict

# Define the GPS epoch
const GPS_EPOCH = DateTime(1980, 1, 6) # DateTime
dt2gpsns(dt) = Dates.value(Millisecond( dt-GPS_EPOCH )) * 1_000_000 # -> Integer
GPS_MORETHAN = dt2gpsns(DateTime(2024,1,1))

# VectorNav data line:
# Tue Mar 28 13:17:11 2023 +0000, 9431257000, -163.267792, -1.067281, -1.441135, -0.011044, 0.011087, -0.989255, 0.145364, 0.000000, 0.000000, 0.000000, 0.234963, -0.005458, 0.101274, 0.432598, -0.534912, -9.797964, 0.000000, 0.000000, 0.000000
#
# 21 science data fields delimited by commas

function read_vecnav_data( file::AbstractString )
    # allocate arrays
    nlmax = countlines(file)
    dt = Vector{DateTime}(undef, nlmax)
    intdata = Matrix{Int64}(undef, nlmax, 2)
    data = Matrix{Float32}(undef, nlmax, 19)

    # read it     >>BLACK PUMAS<<
    nl = 0
    df = dateformat"uuuHH:MM:SSyyyy"
    open( file ) do f
        for line in readlines(f)
            # Tue Mar 28 13:17:11 2023 +0000, 9431257000, -163.267792, -1.067281, -1.441135, -0.011044, 0.011087, -0.989255, 0.145364, 0.000000, 0.000000, 0.000000, 0.234963, -0.005458, 0.101274, 0.432598, -0.534912, -9.797964, 0.000000, 0.000000, 0.000000
            nl += 1
            spl = split(line, r"[\ ,]+")
            try
                dt[nl] = ( 
                    DateTime(spl[2], dateformat"u") + Day(parse(Int8, spl[3])) 
                  + Year(parse(Int16, spl[5])-1) 
                  +(DateTime(spl[4], dateformat"HH:MM:SS") - DateTime(1,1,1)) )
                intdata[nl,:] .= parse.(Int64, spl[6:7])  #  2
                data[nl,:] .= parse.(Float32, spl[8:end]) # 19
            catch
                @show line
            end
        end
    end

    # trim uninitialized times off the start
    # i0 = findfirst(diff(intdata[1:nl,2]) .> 1e17) + 1 # broken
    i0 = findfirst(intdata[1:nl,2] .> GPS_MORETHAN)
    return dt[i0:nl], intdata[i0:nl,:], data[i0:nl,:]
end

function read_vecnav_data( files::AbstractVector )
    # allocate arrays
    nlmax = sum( countlines.(files) )
    dt = Vector{DateTime}(undef, nlmax)
    intdata = Matrix{Int64}(undef, nlmax, 2)
    data = Matrix{Float32}(undef, nlmax, 19)

    # read files
    nl = 0
    df = dateformat"uuuHH:MM:SSyyyy"
    for file in files
        open( file ) do f
            for line in readlines(f)
                # Tue Mar 28 13:17:11 2023 +0000, 9431257000, -163.267792, -1.067281, -1.441135, -0.011044, 0.011087, -0.989255, 0.145364, 0.000000, 0.000000, 0.000000, 0.234963, -0.005458, 0.101274, 0.432598, -0.534912, -9.797964, 0.000000, 0.000000, 0.000000
                spl = split(line, r"[\ ,]+")
                try
                    id = parse.(Int64, spl[6:7])
                    if id[2] > GPS_MORETHAN  # use only initialized GPS times
                        nl += 1 # get data
                        dt[nl] = ( 
                            DateTime(spl[2], dateformat"u") + Day(parse(Int8, spl[3])) 
                        + Year(parse(Int16, spl[5])-1) 
                        +(DateTime(spl[4], dateformat"HH:MM:SS") - DateTime(1,1,1)) )
                        intdata[nl,:] .= id                       #  2
                        data[nl,:] .= parse.(Float32, spl[8:end]) # 19
                    end
                catch
                    @show line
                end
            end
        end
    end

    return dt[1:nl], intdata[1:nl,:], data[1:nl,:]
end

function read_vecnav_dict( file )

    D = Dict{Symbol, Any}()
    vnkeys = Symbol.( [
    "Yaw", "Pitch", "Roll",
    "Quat0", "Quat1", "Quat2", "Quat3",
    "Latitude", "Longitude", "Altitude",
    "MagNED0", "MagNED1", "MagNED2",
    "LinAcc0", "LinAcc1", "LinAcc2",
    "VelNED0", "VelNED1", "VelNED2" ] )

    dt, gpstime, data = read_vecnav_data( file )

    D[:time] = dt
    D[:GpsTime] = gpstime[:,2]
    for (i,k) in enumerate(vnkeys)
        D[k] = data[:,i]
    end
    return D
end

"Read VectorNav data previously concatenated and saved as a JLD2 file. Much faster than reading the text files."
# read_vecnav_dict() = Dict(Symbol(key) => value for (key, value) in load("./data/table/ASTraL_lidarVectorNav.jld2")) # Dict{Symbol, Any}
read_vecnav_dict() = Dict( Symbol(key) => value for (key, value) in
                           load("./data/table/ASTraL_lidarVectorNav.jld2") )

end # module read_vecnav


module chunks

using Dates
using Statistics

export read_stare_time, read_stare_chunk, dt_to_chunkind
export fit_offset

"result is x; set to missing iff i<thr"
masklowi(x, i, thr=1.03) = i<thr ? missing : x

"mean along dimension dims, skipping missing"
missmean(X; dims=1) = mapslices(x -> mean(skipmissing(x)), X, dims=dims)

"anomaly"
anom(x; dims=1) = x.-mean(x; dims=dims)

"find indices i such that each xl[i] is the first >= xs."
function findindices(xs, xl)
    # xs needles define quarries in haystack xl
    xs = filter(x -> x<=last(xl), xs) # prefilter to avoid running off the end of xl
    ind = zeros(Int64, size(xs))
    i = 1
    for (j,x) in enumerate(xs)
        while xl[i] < x
            i += 1
        end
        ind[j] = i
    end
    return ind
end

"average within +-half points of the index of xl"
function indavg(xl, ind; half=5)
    xm = zeros(Float64, size(ind))
    for (i,idx) in enumerate(ind)
        ii = max(1,idx-half) : min(length(xl),idx+half)
        xm[i] = sum(Float64.(xl[ii])) / (2*half+1)
    end
    return xm
end

"""
fit_offset(dt_lidar_chunk) returns offset (seconds) that lidar needs to advance to match the VN,
using hardcoded linear fit derived from pitch solid-body rotation.
"""
function fit_offset(dt)
    slope = -2.4805283261947786e-5
    # 2 pieces
    if dt < DateTime(2024,6,6,11,15)
        origin_time   = DateTime(2024,5,1,23,51,47,117)
        origin_offset = 143.8 #.8062
    else
        origin_time   = DateTime(2024,6,6,23,33,44,321)
        origin_offset = 126.94
    end
    origin_offset + Dates.value(Millisecond(dt - origin_time))/1_000 * slope
end

"""
get_time_shift(mdv, heave) positive result means mdv clock is fast.
sync by subtracting this lag from stare_dt.
"""
function get_time_shift(x, y)
    # filter to make xcorr work better
    xc = xcorr(hp(x[:]), hp(y[:]))
    # plot(-(length(mdv)-1):length(mdv)-1, xc )
    return argmax(xc) - length(x)
end

"return indices of contiguous chunks (separated by < thr) from a series of datetimes"
function dt_to_chunkind(dt, thr=Second(30))
    pickets = findall( t -> t>thr, diff(dt) )
    st_chunk = pickets[1:end-1] .+ 1
    en_chunk = pickets[2:end]
    return st_chunk, en_chunk
end

"hourly files -> chunk time indices"
function read_stare_time( St, dt=Date(2024,5,1) )
    # dt doesn't matter for dividing into chunk indices, as long as St[:time] is contiguous.
    lidar_stare_dt = @. ( DateTime(Date(dt)) + Millisecond(round(Int64, St[:time] * 3_600_000 )) )
    st_chunk, en_chunk = dt_to_chunkind(lidar_stare_dt) # individual stare chunk indices
    # subdivide into shorter chunks???

    # Lidar clock is fast (ahead) by ~60-150 s compared to the GPS, 
    # but don't need to adjust to compute indices of chunks.
    # lidar_clock_fast_by = Millisecond( 1_000 * fit_offset(lidar_stare_dt[st_chunk]) ) # first adjustment
    return st_chunk, en_chunk
end

# "subdivide interval [st en] into fac even intervals"
# subdivide(st,en, fac) = round(Integer, st .+ (st-en)/fac .* [0:fac])
"subdivide single interval [st en] into fac even intervals"
subdivide(st,en, fac) = @. round(Integer, st + ((en-st)*(0:fac)/fac))
# # test
# subdivide(1,240, 4)

# functions to fine sync the 20 Hz VN to the 1 Hz lidar
"true for data that results in finite arithmetic"
good(x) = !ismissing(x) && !isnan(x)

"cov ignoring missing and NaNs"
function goodcov(x, y)
    nn = good.(x) .& good.(y)
    xc = cov(x[nn], y[nn])
end

"covariance between x and yhi[ind] around central indices ind, optionally lagged by r"
xcov(x, yhi, ind=eachindex(yhi), r=0) = goodcov(x, yhi[ind .+ r])
"for iterable rr"
xcov(x, yhi, ind, rr::AbstractVector) = [ xcov(x, yhi, ind, r) for r in rr ] # iterate over offsets

"find optimal offset r within offset_range for fine yhi[ind+r] to sync with x"
function finelagcov(x, yhi, ind, offset_range)
    xc = xcov(x, yhi, ind, offset_range)
    maxcov, fine_offset = findmax(xc)
end

"read and interpolate data to stare chunks"
function read_stare_chunk( dt::TimeType, St, Vn, UV, st, en, ntop=80 )
    # time: truncate the file's datestamp to Date, add the decimal hour
    stare_dt_raw = @. DateTime(Date(dt)) + Millisecond(round(Int64, St[:time] * 3_600_000 )) # 3202
    # synchronize clocks using the previously calculated offset function
    # lidar_clock_fast_by = Millisecond( 126832 ) # adjust for lidar clock fast (ahead) by 126832 milliseconds compared to the GPS.
    lidar_clock_fast_by = Millisecond( round(Int64, 1_000 * fit_offset(stare_dt_raw[st])) )
    stare_dt = stare_dt_raw .- lidar_clock_fast_by # synced to within 1 s
    stare1dt = stare_dt[st:en] # subset

    # pre-subset
    vn_ind = findall( stare_dt[st]-Second(2) .<= Vn[:vndt] .<= stare_dt[en]+Second(2) )

    # fine sync the 20 Hz VN to the 1 Hz lidar
    # St[:pitch][st:en] #  1 Hz
    # Vn[:roll][vn_ind] # 20 Hz
    ind0 = findindices( stare1dt, Vn[:vndt][vn_ind] ) # main indices for VN
    # fine_offset_range = -20:20 # offset indices
    # maxcov, fine_offset = finelagcov(St[:pitch][st:en], Vn[:Roll][vn_ind], ind0, fine_offset_range)
    # fine_offset seems to introduce phase errors
    
    # dopplervel (time, z) masked by intensity threshold
    dopplervel = masklowi.(St[:dopplervel][st:en,1:ntop], St[:intensity][st:en,1:ntop])
    # mdv = missmean(dopplervel, dims=2)[:] # conditional mean can have biases
    jj = all( isfinite.(dopplervel), dims=1) # heights with all good data
    mdv = mean(dopplervel[:,jj], dims=2)[:] # mean of filled heights

    # interpolate Ur,Vr, heave to the lidar stare grid
    # ind = ind0 .+ fine_offset
    ind = ind0 #.- fine_offset
    pitch = indavg( Vn[:Pitch][vn_ind], ind) # 11-point centered mean around ind
    roll  = indavg( Vn[:Roll ][vn_ind], ind)
    VelNED0 = indavg( Vn[:VelNED1][vn_ind], ind) # note VN mounted 90-degrees off
    VelNED1 = indavg( Vn[:VelNED0][vn_ind], ind)
    VelNED2 = indavg( Vn[:VelNED2][vn_ind], ind) # +down
    # # resync the clock to the VactorNav heave - brittle
    # ind = findindices( stare1dt, Vn[:vndt] )
    # heave = Vn[:VelNED2][ind] # nearest neighbor interp to lidar times
    # shift = get_time_shift(mdv[:],heave[:]) # to nearest second, superfluous.
    # # interpolate for the updated synced time
    # stare1dt .-= Millisecond((1020-80)*shift)
    # ind = findindices( stare1dt, Vn[:vndt] ) # this works
    # heave = indavg( Vn[:VelNED2], ind)

    # mean relative velocity
    Ur = zeros(size(dopplervel))
    Vr = zeros(size(dopplervel))
    ind = findindices( Dates.value.(stare1dt), Dates.value.(UV["time"]))
    # result must be 1:1 for stare1dt and ind
    ls = length(stare1dt)
    li = length(ind)
    if li < ls # extend ind with last index of UV
        ind = [ind; length(UV["time"]).+zeros(Int32, ls-li)]
    end
    for ih in 1:ntop # loop to broadcast to consistent size
        Ur[:,ih] .= UV[:ur][ind,ih]
        Vr[:,ih] .= UV[:vr][ind,ih]
    end

    # questionable: fill all the mean relative velocities
    isf = isfinite.(Vr)
    Vr[.!isf] .= mean(Vr[isf])
    Ur[.!isf] .= mean(Ur[isf])
    
    return dopplervel, pitch, roll, VelNED0, VelNED1, VelNED2, Ur, Vr, mdv
end

end # module chunks