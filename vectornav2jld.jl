## load the VectorNav data and save as a JLD2 file

using Revise
using Pkg; Pkg.activate(".")

using Dates
# using Statistics
# using StatsBase
using Interpolations
using JLD2

include("./readers.jl")
using .ShipPosmv, .NoaaDas

include("./read_lidar.jl")
using .read_vecnav: read_vecnav_dict

include("./timing_lidar.jl")
using .timing_lidar

## helper functions
pd = permutedims
m2n(x) = ismissing(x) ? NaN : x

anom(x; dims=1) = x .- mean(x; dims=dims)

# functional wrapper for PyPlot plot
plotf(f::Function, x, va...) = plot(x, f.(x), va...)

"""
Concatenate Array variables from two Dicts into one
along the dimension dim
that are values of the same key of dict1 and dict2.
"""
function cat_dicts(dict1::Dict, dict2::Dict; dim=1)
    result_dict = typeof(dict1)()
    for key in intersect( keys(dict1), keys(dict2) )
        result_dict[key] = cat(dict1[key], dict2[key], dims=dim)
    end
    return result_dict
end

## precompile reader methods on a short file
test = read_vecnav_dict( joinpath(homedir(), "Data/EKAMSAT/lidar/table/VNshort1.txt") )
test = read_vecnav_dict( joinpath.(homedir(), ["Data/EKAMSAT/lidar/table/VNshort1.txt", "Data/EKAMSAT/lidar/table/VNshort1.txt"]) )

## read VectorNav data
# leg 1
datadir = joinpath(homedir(), "Data/EKAMSAT/lidar/table/")
file = joinpath(datadir, "leg1", "VectorNavTableData.txt")

@time Vn1 = read_vecnav_dict( file ) # 97 s on Mac ARM laptop, 147 s on Windoz Intel

# leg 2
file = joinpath(datadir, "leg2", "VectorNavData2023_04_30.txt") # ignore wrong datestamp on filename
@time Vn2 = read_vecnav_dict( file ) # 65 s
# shorter than the leg 1 file because only contains 2nd half of leg 2

"convenience wrapper for making a nice time axis from time and gpstime"
f(D::Dict) = timing_lidar.gpstime2gpsvndt( D[:time], D[:GpsTime] )
vndt1 = f( Vn1 )
vndt2 = f( Vn2 )
vndt = [vndt1; vndt2]

## concatenate data from the 2 leg files
Vn = cat_dicts(Vn1, Vn2)
Vn[:vndt] = vndt
# save as a JLD2 file for faster recall
save("./data/table/ASTraL_lidarVectorNav.jld2", Dict(String(key) => value for (key, value) in Vn))

# load Vn Dict, takes 5 s, converts keys back to Symbols
Vn = Dict(Symbol(key) => value for (key, value) in load("./data/table/ASTraL_lidarVectorNav.jld2"))
