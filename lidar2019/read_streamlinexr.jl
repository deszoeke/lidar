## preamble

using Pkg; Pkg.activate(".")

using Revise
using Dates
# using CSV
# using DataFrames
using DelimitedFiles
using Printf


## functions for reading hpl files

"""
header = read_streamlinexr_head(file_path)
Read the Halo Photonics StreamLineXR .hpl file header.
"""
function read_streamlinexr_head(file_path)
    # put the header in a Dict
    h = Dict(
        :nlines => countlines(file_path),
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
modifying read_streamlinexr_data!(file_path, header, beams)
Read data and fill in the beams.
"""
function read_streamlinexr_data!(file_path, h, beams)
    # use header information in h
    nlines = h[:nlines]
    ngates = h[:ngates]

    # beams could be rays or times
    nbeams = round(Int, (nlines-17) / (1+ngates)) # = nrays*ntimes
    beam_timeangles = zeros(nbeams, 5)
    beam_velrad = zeros(nbeams, ngates, 4)

    # for User wind profiles beam <--> VAD ray
    # for Stare beam <--> time

    # open and read the file
    open(file_path) do file
        for _ in 1:17 # skip 17 header lines
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
    beams[:time     ] .= beam_timeangles[:,1] # decimal hours
    beams[:azimuth  ] .= beam_timeangles[:,2] # degrees
    beams[:elevangle] .= beam_timeangles[:,3] # degrees
    beams[:pitch    ] .= beam_timeangles[:,4]
    beams[:roll     ] .= beam_timeangles[:,5]
    # by gate
    beams[:height   ] .= (beam_velrad[1,:,1].+0.5) .* h[:gatelength] # center of gate, assumes same for all beams

    # dependent variables (beam, gate)
    beams[:dopplervel] .= beam_velrad[:,:,2] # m/s
    beams[:intensity ] .= beam_velrad[:,:,3] # SNR + 1
    beams[:beta      ] .= beam_velrad[:,:,4] # m-1 sr-1  backscatter?
end

"""
beams, h = read_streamlinexr(file_path)
Read the Photonics StreamLineXR .hpl file data and header.
"""
function read_streamlinexr(file_path)
    # use header information in h
    h = read_streamlinexr_head(file_path)
    nlines = h[:nlines]
    ngates = h[:ngates]

    # beams could be rays or times
    nbeams = round(Int, (nlines-17) / (1+ngates)) # = nrays*ntimes
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
    read_streamlinexr_data!(file_path, h, beams)
    return beams, h
end


## test read one file
file_path = "./data/Stare_06_20190710_06.hpl"
beams, header = read_streamlinexr(file_path)