#using MAT
using NCDatasets
using Dates
using Interpolations

using PyPlot
#using Plots, UnicodePlots

using PyCall
using PyCall: PyObject

# allow for plotting with missing values
function PyObject(a::Array{Union{T,Missing},N}) where {T,N}
    numpy_ma = PyCall.pyimport("numpy").ma
    pycall(numpy_ma.array, Any, coalesce.(a,zero(T)), mask=ismissing.(a))
end

pd = permutedims

#file = matopen("../2020-01-04/microDop_2020-01-04.10.mat")
#data = read(file, "mdData")["data"] # returns Dict

# lidar data
datapaths = ["/mnt/c/Users/deszoeks/Data/", joinpath(homedir(),"Data"), joinpath(homedir(),"data")]
datapath = filter(isdir, datapaths)[1]
velpath = joinpath(datapath, "ATOMIC/microDop/reprocessed20210224/")
wndpath = joinpath(datapath, "ATOMIC/microDop/windprofiles/")
navpath = joinpath(datapath, "ATOMIC/flux")

#prototypefile = "EUREC4A_ATOMIC_RonBrown_lidar_microDop02_lev1_20200113_0100-20200113_0200.nc"
prefix = "EUREC4A_ATOMIC_RonBrown_lidar_microDop02_lev1_"
extn = ".nc"
datafilelist = readdir(joinpath(datapath,velpath)) |> filter(startswith(prefix)) |> filter(endswith(extn))

# parse dates in filenames
todt(d::T where T<:AbstractString,t::S where S<:AbstractString) = DateTime(d*t, dateformat"yyyymmddHHMM")
#todt(d::S,t::P) = DateTime(d*t, dateformat"yyyymmddHHMM") where S<:AbstractString where P<:AbstractString
filestartdate(filename) = todt(split( filename, ('_','-','.') )[end-4:end-3]...)
fileenddate(filename)   = todt(split( filename, ('_','-','.') )[end-2:end-1]...)
dtstarts = filestartdate.(datafilelist)
dtends   = fileenddate.(  datafilelist)

getfilename(dt) = joinpath( datapath, velpath, datafilelist[dtstarts .<= dt .< dtends][1] )
lidardataset(dt) = NCDataset( getfilename(dt) )

ds = lidardataset(DateTime(2020,1,7,11,30))

yearsec_to_dt(ys, year=2020) = DateTime(year) + Millisecond(round(Integer, 1e3*ys))

range = ds[:scannerHeight] .+ filter(isfinite, ds[:range][:])
nrange = length(range)

#cond(snr) = !ismissing(snr) && isfinite(snr) && snr >= -26
#m = cond.(ds[:snrData][1:nrange,:])
#maskit(x, m=m) = m ? x : missing
nanmask(x, m) = m ? x : NaN # masking with NaNs works for PyPlot
ys = ds[:time][:]
w   = ds[:mcVelData][1:nrange,:]
snr = ds[:snrData][1:nrange,:]

# plot
#=
clf()
subplot(2,1,1)
pcolormesh( yearsec_to_dt.(ys), range./1e3, nanmask.(w,snr .>= -25 ),
	   vmin=-1, vmax=1, cmap=ColorMap("RdBu_r") )
colorbar()
ylim([0, 2.5])
PyPlot.title("Doppler velocity (m/s)")
PyPlot.ylabel("height (km)")

subplot(2,1,2)
pcolormesh( yearsec_to_dt.(ys), range./1e3, snr, vmax=5, vmin=-26, cmap=ColorMap("RdBu_r") )
colorbar()
ylim([0, 2.5])
PyPlot.title("signal to noise ratio (dB)")
PyPlot.ylabel("height (km)")
PyPlot.xlabel("time")

tight_layout()
PyPlot.savefig("vel_snr.png")
#PyPlot.savefig("vel_snr.svg")
=#

#UnicodePlots(heatmap(pd(data["velocity"][1:1000][1:200]))

# mean wind data
#wndpath = "/mnt/c/Users/deszoeks/Data/ATOMIC/microDop/windprofiles"
windfile = "EUREC4A_ATOMIC_RonBrown_lidar_microDop02_lev2_windProfs_RTR_20200107_2100-20200211_2300.nc"
wnds = NCDataset( joinpath(wndpath, windfile) )

meanys = wnds[:time][:]
tt = isfinite.(meanys)
meanheight = wnds[:height][:]
hh = isfinite.(meanheight)

meanwnddt = yearsec_to_dt.(meanys)
# winds out of NE: U,V both negative
U = @. -wnds[:speed][:,:] * sind(wnds[:direction][:,:]) # + to east
V = @. -wnds[:speed][:,:] * cosd(wnds[:direction][:,:]) # + to north

# interpolate using second of year
mint(U) = interpolate((meanheight[hh], meanys[tt]), U[hh,tt], Gridded(Linear()))
mintU = mint(U) # explicitly define interpolator methods before calculating!
mintV = mint(V)
# interpolate to middle of chunks
chmid = @. Dates.value(Second((DateTime(2020,1,9,0,2,30):Minute(5):meanwnddt[end]) - DateTime(2020,1,1,0)))

mU = zeros(100, length(chmid))
mV = zeros(100, length(chmid))
for ci in CartesianIndices(mU)
    mU[ci] = mintU(range[ci[1]], chmid[ci[2]])
    mV[ci] = mintV(range[ci[1]], chmid[ci[2]])
end

# need RELATIVE wind!
# read ship navigaton data
navfile = "EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc"
navds = NCDataset( joinpath(navpath, navfile) )
navys = navds[:time][:] # second of year
Base.deg2rad(x::Missing) = missing # method for missing data
unav = @. navds[:sog][:] * sind.(navds[:cog][:])
vnav = @. navds[:sog][:] * cosd.(navds[:cog][:])
isdata(x) = !ismissing(x) && isfinite(x)
nn = @. isdata(unav) && isdata(vnav) && isdata(navys) && navds[:ship_contamination_flag]==0
unavint = interpolate((collect(skipmissing(navys[nn])),), collect(skipmissing(unav[nn])), Gridded(Linear())) # interpolation method
vnavint = interpolate((collect(skipmissing(navys[nn])),), collect(skipmissing(vnav[nn])), Gridded(Linear()))
Unavi = unavint.(chmid) # calc the interpolation
Vnavi = vnavint.(chmid)

Urel = mU .- pd(Unavi)
Vrel = mV .- pd(Vnavi)
Srel = @. sqrt( Urel^2 + Vrel^2 )

# this how the relative wind variable is arranged
mdates = PyPlot.matplotlib.dates
clf()
fig, ax = subplots(1,1)
pcm = ax.pcolormesh(DateTime(2020,1,1).+Second.(chmid), range[1:100]/1e3, Srel); colorbar(pcm)
ax.xaxis.set_major_formatter( mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()) )
ax.set_ylabel("height (km)")
ax.set_title("ship-relative horizontal wind\ninterpolated to center of 5-min windows")
fig.tight_layout()
