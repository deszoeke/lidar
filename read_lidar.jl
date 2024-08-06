# using Revise
# using Pkg; Pkg.activate("/Users/deszoeks/Projects/ASTRAL/lidar")

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

### utilites
pd = permutedims
m2n(x) = ismissing(x) ? NaN : x

### read stares functions
# functions for reading hpl files
lidardir = "./data/lidar"

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
function read_streamlinexr_stare!(file_path, h, beams, nheaderlines=17; nbeams0=0)
    # use header information in h
    nlines = h[:nlines]
    ngates = h[:ngates]

    # beams could be rays or times
    nbeams = round(Int, (nlines-nheaderlines) / (1+ngates)) # = nrays*ntimes
    beam_timeangles = zeros(nbeams, 5)
    beam_velrad = zeros(nbeams, ngates, 4)

    # for User wind profiles beam <--> VAD ray
    # for Stare beam <--> time

    # open and read the file
    open(file_path) do file
        for _ in 1:nheaderlines # skip header lines
            readline(file)
        end

        # now read data
        for ibeam = 1:nbeams
            # beam described by a batch of 1+ngates lines
            # Read the beam parameter line
            line = readline(file)
            beam_timeangles[ibeam,:] .= parse.(Float64, split(line))
            # Read each gate in the beam
            for igate = 1:ngates
                line = readline(file)
                beam_velrad[ibeam, igate,:] .= parse.(Float64, split(line))
            end
        end
    end # close the file

    # parse the variables into the dict beams
    # by beam
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

function read_streamlinexr_stare(file_path::AbstractString, nheaderlines=17)
    # Not implemented: could loop over a vector of files.

    # use header information in h 
    h = read_streamlinexr_head(file_path)
    nlines = h[:nlines]
    ngates = h[:ngates]

    # beams could be rays or times
    nbeams = round(Int, (nlines-nheaderlines) / (1+ngates)) # = nrays*ntimes
    # initialize a beams Dict
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
    read_streamlinexr_stare!(file_path, h, beams, nheaderlines)
    return beams, h
end

"read multiple files"
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
        )

    # read each file and fill beams with data
    nbeams0 = 0
    h = read_streamlinexr_head(file_path[1]) # first header
    for (i,file) in enumerate(file_path)
        h = read_streamlinexr_head(file) # reread header for each file
        read_streamlinexr_beam_timeangles!(file, h, beams, nheaderlines; nbeams0=nbeams0) # updates beams[keys][nbeams0 .+ (1:nbeams[i])]
        nbeams0 += nbeams[i] # number of beams now read
    end
    return beams, h, nbeams0
end

### read mean wind functions
uvdir = "./data/lidar/netcdf/"

"get daily mean-wind netcdf dataset"
function get_daily_meanuv( dt::Union{Date, DateTime} )
    get_daily_meanuv( Dates.format(Date(dt), dateformat"yyyymmdd") )
end
function get_daily_meanuv( yyyymmdd::AbstractString )
    Dataset(joinpath(uvdir,"ekamsat_lidar_uv_$(yyyymmdd).nc"))
end

### read VectorNav functions
Vndir = "./data/lidar/table/" # uses symbloic link ./data in cwd

"read daily JLD2 as a Dict"
function read_daily_Vn( dt::Union{Date, DateTime} )
    read_daily_Vn( Dates.format(Date(dt), dateformat"yyyymmdd") )
end
function read_daily_Vn( yyyymmdd::AbstractString )
    Vn = JLD2.load(joinpath(Vndir, "VectorNavTable_$(yyyymmdd).jld2"))
end

end # module read_lidar