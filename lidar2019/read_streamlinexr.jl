using CSV
using DataFrames

function read_streamlinexr(file::String, nrange::Int)
    # Read data from file skipping header
    x = CSV.read(file, DataFrames.header=false, DataFrame, limit=10000, footerskip=4)

    # get ray time and angles
    # Decimal time (hours)  Azimuth (degrees)  Elevation (degrees) Pitch (degrees) Roll (degrees)
    timeangles = x[1:nrange:end, :]
    a = Dict(
        :dechour => timeangles[!, 1],
        :azimuth => timeangles[!, 2],
        :elevation => timeangles[!, 3],
        :pitch => timeangles[!, 4],
        :roll => timeangles[!, 5]
    )

    nt = length(a[:dechour])

    # get ray data
    idata = mod.(0:size(x, 1)-1, nrange+1) .> 0
    a[:rangegate ] = reshape(x[idata, 1], (nrange, nt))
    a[:dopplervel] = reshape(x[idata, 2], (nrange, nt))
    a[:intensity ] = reshape(x[idata, 3], (nrange, nt))
    a[:beta      ] = reshape(x[idata, 4], (nrange, nt))

    return a
end
