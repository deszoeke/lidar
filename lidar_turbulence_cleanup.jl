# preamble
using Revise
using Pkg; Pkg.activate(".")

using Dates
using Statistics
using Rotations
using Interpolations
using DSP
using FFTW
using NCDatasets
using JLD2
using Printf

include("./read_lidar.jl")
using .read_lidar
using .read_lidar.stare
using .read_vecnav: read_vecnav_dict
import .chunks
# explicitly load into Main global scope
read_stare_time  = Main.chunks.read_stare_time
# read_stare_chunk = Main.chunks.read_stare_chunk # redefined below
fit_offset = Main.chunks.fit_offset

include("./timing_lidar.jl")
using .timing_lidar
include("./readers.jl")
using .NoaaDas: cat_dicts
# using MAT

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

# redefine read_stare_chunk
"""
Motion from Doppler velocity, by inferring displacement of the median velocity 
by the platform motion at each time.
"""
function get_mdv(f, intensity, dopplervel)
    # levels with Doppler data at all times in the chunk
    goodlevels = findall(all(f, intensity; dims=1)[:]) # bool --> indices
    # median over whole chunk, continuously available levels
    chunkmedian = median(dopplervel[:,goodlevels]) # over both dims
    # median over vertical dim at each time
    timeslicemedian = median(dopplervel[:,goodlevels], dims=2)
    mdv = timeslicemedian .- chunkmedian
end

"read and interpolate data to stare chunks"
function read_stare_chunk( dt::TimeType, St, Vn, UV, st, en, ntop=80 )
    noisethr = 1.03
    f(x) = !ismissing(x) && x > noisethr

    # 2025 09 20
    stare_dt_raw = @. DateTime(Date(dt)) + Millisecond(round(Int64, St[:time][st:en] * 3_600_000 )) # 3202
    # synchronize clocks using the previously calculated offset function
    # lidar_clock_fast_by = Millisecond( 126832 ) # adjust for lidar clock fast (ahead) by 126832 milliseconds compared to the GPS.
    lidar_clock_fast_by = Millisecond( round(Int64, 1_000 * fit_offset(stare_dt_raw[1])) )
    stare_dt = stare_dt_raw .- lidar_clock_fast_by # synced to within 1 s
    stare1dt = stare_dt[:] # subset

    # dopplervel (time, z) masked by intensity threshold
    dopplervel = masklowi.(St[:dopplervel][st:en,1:ntop], St[:intensity][st:en,1:ntop])
    pr_(v) = ismissing(v) ? v : v[:]
    mdv = pr_( get_mdv(f, St[:intensity][st:en,1:ntop], dopplervel) )

    # pre-subset
    vn_ind = findall( stare1dt[1]-Second(2) .<= Vn[:vndt] .<= stare1dt[end]+Second(2) )
    if length(vn_ind) > 0
        ind = findindices( stare1dt, Vn[:vndt][vn_ind] ) # main indices for VN
        pitch = indavg( Vn[:Pitch][vn_ind], ind) # 11-point centered mean around ind
        roll  = indavg( Vn[:Roll ][vn_ind], ind)
        VelNED0 = indavg( Vn[:VelNED1][vn_ind], ind) # note VN mounted 90-degrees off
        VelNED1 = indavg( Vn[:VelNED0][vn_ind], ind)
        VelNED2 = indavg( Vn[:VelNED2][vn_ind], ind) # +down
    else # no motion data
        pitch = roll = VelNED0 = VelNED1 = VelNED2 = fill(missing, length(stare1dt))
    end

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
        Ur[:,ih] .= UV[:ur][ih, ind] # 2025 09 20 fixed swapped subscripts
        Vr[:,ih] .= UV[:vr][ih, ind]
    end

    # questionable: fill all the mean relative velocities
    isf = isfinite.(Vr)
    Vr[.!isf] .= mean(Vr[isf])
    Ur[.!isf] .= mean(Ur[isf])
    
    return dopplervel, pitch, roll, VelNED0, VelNED1, VelNED2, Ur, Vr, mdv
end

# function library with utility functions,  functions for subsetting, for displacements, and for structure functions

# utility functions
pd = permutedims
m2n(x) = ismissing(x) ? NaN : x
good(x) = !ismissing(x) & isfinite(x)

"""
    binavg(y, x, b; f=identity, w=y->1)
Bin average y(x) in bins b of coordinate x.
Skip missing by passing the optional function arguments
f(y) = ismissing(y) ? 0 : y;   w(y) = !ismissing(y)
"""
function binavg(y, x, b; f=identity, w=y->1)
    a = zeros(length(b))
    c = zeros(length(b))
    for (i,x) in enumerate(x)
        bi = findlast(j -> j < x, b) # findlast(<(x), b)
        a[bi] += f(y[i])
        c[bi] += w(y[i])
    end
    return a./c
end

# functions for masking and averaging data

"NaN -> missing"
n2m(x) = isfinite.(x) ? x : missing

"result is x; set to missing iff i<thr"
masklowi(x, i, thr=1.03) = i<thr ? missing : x

"mean along dimension dims, skipping missing"
missmean(X; dims=1) = mapslices(x -> mean(skipmissing(x)), X, dims=dims)

noisethr = 1.03
f(x) = x>noisethr # use intensity for x

# # mean Doppler velocity # OLD and tested
# function get_mdv(f, dopplervel)
#     X = missingnoise.(f, dopplervel)
#     mapslices(x->mean(skipmissing(x)), X; dims=2)[:]
# end

# Suggest only computing mdv (for motion removal) where
# there are nonnoise DVs for the whole chunk, so not to alias sampling differences.
# Then, only take the mean of the middle 5th quantile, to smoothly approximate the mode

"""
Motion from Doppler velocity, by inferring displacement of the median velocity 
by the platform motion at each time.
"""
function get_mdv(f, intensity, dopplervel)
    # levels with Doppler data at all times in the chunk
    goodlevels = findall(all(f, intensity; dims=1)[:]) # bool --> indices
    if length(goodlevels) < 5
        mdv = missing # not enough levels with data
    else
        # median over whole chunk, continuously available levels
        chunkmedian = median(dopplervel[:,goodlevels]) # over both dims
        # median over vertical dim at each time
        timeslicemedian = median(dopplervel[:,goodlevels], dims=2)
        mdv = timeslicemedian .- chunkmedian
    end
    return mdv
end
# UNTESTED!!
# This modification circumvents issues with sampling and better separates the motion estimate
# from the air velocities, assuming the median of the overall velocity distribution is 
# unchanged, and the Doppler velocity is merely shifted at each time by the platform motion,
# but with zero average shift since mean(motion) = median(motion) = 0.

# Remove median Doppler velocity.
# Still plot dopplervel even if mdv scalar missing
# prevents dissipation calculation.
function remove_mdv(f, intensity, dopplervel)
    mdv = get_mdv(f, intensity, dopplervel)
    if ismissing(mdv)
        return dopplervel # not enough good data to compute mdv
    else
        return dopplervel .- mdv
    end
end
function remove_mdv!(f, intensity, dopplervel)
    mdv = get_mdv(f, intensity, dopplervel)
    if !ismissing(mdv)
        dopplervel .-= mdv
    end
end

"anomaly"
anom(x; dims=1) = x.-mean(x; dims=dims)

# highpass filter
"""
hp(x, fcutoff=1/80)    highpass filter x,
by default filtfilt 4th-order Butterworth, fs=1
"""
function hp(x, fcutoff=1/80;
    order=4,
    designmethod=Butterworth(order), 
    fs=1,
    responsetype = Highpass(fcutoff; fs=fs) )
    
    filtfilt(digitalfilter(responsetype, designmethod), x)
end


# make simple linear temporal interpolation
# maybe fast
# most time is spent searching for indices
# indices are monotonic

"""
Return all the indices i such that each xl[i] is the first >= each xs.
Assumes xs, xl are ordered and loops through xs only once.
Quarry for needles xs in haystack xl.
"""
function findindices(xs, xl)
    # needles xs define quarries in haystack xl
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

"average xl within windows to right of points of the index ind of xl"
function indavg(xl, ind; full=20)
    xm = zeros(Float64, size(ind))
    for (i,idx) in enumerate(ind)
        ii = max(1,idx) : min(length(xl),idx+full)
        # xm[i] = sum(Float64.(xl[ii])) / (full+1)
        xm[i] = mean(Float64.(xl[ii]))
    end
    return xm
end

# test data (precompiles)
let xl = 1:60_000_000, xs = 20:20:60_000_000
    ind = findindices(xs, xl)
    indavg(xl, ind)
end

# Adjust true vertical velocity for relative wind * sin(tilt)
# and the platform velocity
trigs(pitch, roll) = ( cos(pitch), sin(pitch), cos(roll), sin(roll) )
# cospitch, sinpitch, cosroll, sinroll = trigs(pitch, roll)

function wtrue_trigs(w, Ur, Vr, pitch, roll)
    # approximate, better to use rotations
    cospitch, sinpitch, cosroll, sinroll = trigs(pitch, roll)
    wtrue = ( w + Ur*sinpitch - Vr*cospitch*sinroll ) / (cospitch*cosroll)
end

"""
wtrue(dopplervel, Ur, Vr, heaveveldown, roll, pitch)
Return true radial velocity component in lidar beam frame (+away).
Rotate vertical VelNED and mean ship-relative wind (Ur, Vr)
from inertial level ship coorindates 
to lidar beam coordinates using roll and pitch.
"""
function wtrue( dopplervel, Ur, Vr, heaveveldown, roll, pitch )
    # external ship frame
    vvn_ship = [0, 0, heaveveldown] # VectorNav vertical velocity vector (NED coordinate)
    wnd_ship = [Ur, Vr, 0]          # mean horizontal relative wind, w=0 (NED coordinate)
    wnd_vn_ship = wnd_ship - vvn_ship

    # rotate from ship NED frame to lidar NED frame
    R = RotX(roll*π/180) * RotY(pitch*π/180)

    # mean vertical-radial-lidar relative velocity in the lidar platform body frame (NED)
    # includes heave-induced velocity
    wnd_lidar =  R * wnd_vn_ship # lidar NED frame (down-positive) vector

    # signs: lidar upward heave vel > 0 ==> lidar VelNED2 < 0, induced radial velocity < 0 (towards)

    # scalar true radial velocity (+up), adjusting for heave velocity
    # and mean wind component in beam direction.
    # wturb and dopplervel is away-positive. true radialvel is dopplervel + platform vel
    # wtrue = wrel + wplatform
    # trueradialvel is +up
    trueradialvel = dopplervel + -wnd_lidar[3] # negate downward wnd_lidar: NED +down, dopplervel +up
end

function wtrue( dopplervel, surgevel, swayvel, heaveveldown, Ur, Vr, roll, pitch )
    # external ship frame
    vvn_ship = [surgevel, swayvel, heaveveldown] # VectorNav vertical velocity vector (NED coordinate)
    #vvn_ship = [VelNED0, VelNED1, VelNED2]
    wnd_ship = [Ur, Vr, 0]          # mean horizontal relative wind, w=0 (NED coordinate)
    wnd_vn_ship = wnd_ship - vvn_ship

    # rotate from ship NED frame to lidar NED frame
    R = RotX(roll*π/180) * RotY(pitch*π/180)

    # mean vertical-radial-lidar relative velocity in the lidar platform body frame (NED)
    # includes heave-induced velocity
    wnd_lidar =  R * wnd_vn_ship # lidar NED frame (down-positive) vector

    # signs: lidar upward heave vel > 0 ==> lidar VelNED2 < 0, induced radial velocity < 0 (towards)

    # scalar true radial velocity (+up), adjusting for heave velocity
    # and mean wind component in beam direction.
    # wturb and dopplervel is away-positive. true radialvel is dopplervel + platform vel
    # wtrue = wrel + wplatform
    # trueradialvel is +up
    trueradialvel = dopplervel + -wnd_lidar[3] # negate downward wnd_lidar: NED +down, dopplervel +up
end

# functions for indexing sample pairs for structure functions

# generate unique pairs of indices
"index unique pairs in a vector of length n"
function uniquepairs(n) 
    [ [l1, l2] for l1 in 1:n for l2 in (l1+1):n ]
end
"index pairs of points in adjacent levels"
allcross(n) = [ [l1, l2] for l1 in 1:n for l2 in 1:n ]

# beam geometry
"lidar beam range"
rng(iz, rangegate=24.0) = rangegate * (iz-1 + 0.5)

"""
compile indices of lidar volumes to be compared with
structure functions
"""
function lidarindices(nt, nz, z1=1; nlevelstats=1)
    if nlevelstats == 3
        # The complete set that doesn't repeat pairs is 
        # 1 the complete set of nt*(n-1)/2 pairs for the top level (3)
        # 2 the 2*nt*nt sets of pairs between every point in top (3) level and the next 2 levels
        # Iteratively slide this box upward by 1 level for each level.
    
        # index pairs in middle level 2-2
        up = uniquepairs(nt)
        it1 = map(i->i[1], up) # time indices for pairs of point1, point2
        it2 = map(i->i[2], up)
        ci1_r22 = CartesianIndex.(tuple.(it1,z1)) # 1st point in pair lev
        ci2_r22 = CartesianIndex.(tuple.(it2,z1)) # 2nd 
    
        # index pairs of points from level 2-1, and 2-3
        ac = allcross(nt)
        it1 = map(i->i[1], ac)
        it2 = map(i->i[2], ac)
        ci1_r21 = ci1_r23 = CartesianIndex.(tuple.(it1,2))
        ci2_r21 = CartesianIndex.(tuple.(it2,z1-1))
        ci2_r23 = CartesianIndex.(tuple.(it2,z1+1))
    
        # omnibus set of cartesian index pairs for a level, including points in lev above and below
        ci1 = [ci1_r23; ci1_r22; ci1_r21] # first of pairs
        ci2 = [ci2_r23; ci2_r22; ci2_r21]
        li1 = LinearIndices(ci1)
        li2 = LinearIndices(ci2)
        
    elseif nlevelstats == 1
        # just use structure function velocity pairs from one level of lidar range
        up = uniquepairs(nt)
        it1 = map(i->i[1], up) # time indices for pairs of point1, point2
        it2 = map(i->i[2], up)
        ci1_r11 = CartesianIndex.(tuple.(it1,z1)) # 1st point in pair lev
        ci2_r11 = CartesianIndex.(tuple.(it2,z1)) # 2nd point in same lev
    
        # set of cartesian index pairs for a level, including points in lev above and below
        ci1 = ci1_r11 # first of pairs
        ci2 = ci2_r11
        li1 = LinearIndices(ci1)
        li2 = LinearIndices(ci2)
    end
    
    it1 = map(idx->idx[1], ci1) #  t index of first point(s)
    iz1 = map(idx->idx[2], ci1) #  z index of first
    it2 = map(idx->idx[1], ci2) #  t       of second points(s)
    iz2 = map(idx->idx[2], ci2) #  z          second
    
    return ci1,ci2, li1,li2, it1,iz1,it2,iz2
end

# try example
ci1,ci2, li1,li2, it1,iz1,it2,iz2 = lidarindices(1000, 80)

# functions for displacments and structure functions 

rangegate = 24.0 # m # for ASTRAL 2024 Halo Photonics StreamLineXR
timestep = 1.02 # s # time between samples

"""
zm, dr2, dz2, D2 = displacements( ci1,ci2, Udt,Vdt, pitch,roll, w; rangegate=rangegate)
Displacements of sample pairs for one (vertical) subvolume.
"""
function displacements( ci1,ci2, Udt,Vdt, pitch,roll, w; rangegate=rangegate , timestep=timestep)
    # get the individual indices
    it1 = map(idx->idx[1], ci1) #  t index of first point(s)
    iz1 = map(idx->idx[2], ci1) #  z index of first
    it2 = map(idx->idx[1], ci2) #  t       of second points(s)
    iz2 = map(idx->idx[2], ci2) #  z          second

    rng(iz) = rangegate * (iz-1 + 0.5) # center of gates

    # horiz translation of the sample volumes by mean wind
    Udtbar = @. (Udt[iz2] + Udt[iz1]) / 2
    Vdtbar = @. (Vdt[iz2] + Vdt[iz1]) / 2
    X = @. Udtbar * (it2 - it1)
    Y = @. Vdtbar * (it2 - it1)
    # vertical middle of pair
    zm = @. (rng(iz2) * cos(pitch[it2])*cos(roll[it2]) + rng(iz1) * cos(pitch[it1])*cos(roll[it1])) / 2
    # displacement between pair of points
    dz = @.     rng(iz2) * cos(pitch[it2])*cos(roll[it2]) - rng(iz1) * cos(pitch[it1])*cos(roll[it1])
    dx = @. X + rng(iz2) *-sin(pitch[it2])                - rng(iz1) *-sin(pitch[it1])
    dy = @. Y + rng(iz2) * cos(pitch[it2])*sin(roll[it2]) - rng(iz1) * cos(pitch[it1])*sin(roll[it1])
    # distance between
    dz2 = dz .* dz
    dr2 = @. dz2 + dx*dx + dy*dy
    # CORRECT W for HEAVE and for TILTING into the horizontal wind
    # vel structure function
    D2 = @. (w[ci2] - w[ci1])^2
    # return properties of pairs
    return zm, dr2, dz2, D2
end

"dr^2/3 (1-(dz/dr)^2/4) displacement function for computing dissipation from structure function pairs"
rhopair(dr2, dz2) = dr2^(1/3) * (1 - dz2/(4*dr2))

# structure function dissipation functions

# stucture function constants
C2ll = 2.0
epsilon(A) = sqrt(3/4 * A/C2ll)^3
# struf(epsilon, r,r1) = C2ll * epsilon^(2/3) * r^(2/3) * (4 - (r1/r)^2)/3
# instruf(w1,w2) = (w1-w2)^2
# rho(r1,r) = r^(2/3) * (1 - ((r1/r)^2)/4)
# zmid(z1,z2) = (z1 + z2) / 2
# plot bin averaged instruf vs rho
# fit 
# D = A*rho + noise
# for A and noise
# A = 4/3 * C2ll * epsilon^(2/3)

"bin average D2 in equally-populated bins of rho"
function equal_bin(rho, D2; nbin=200, nbin_out_max=17 )
    ii = findall(.!ismissing.(rho) .& .!ismissing.(D2) )
    nrho = length(ii)
    if nrho >= 20
        sp = sortperm(rho[ii])
        srho = rho[ii][sp]
        step = max(1,round(Int32,nrho/nbin))
        rhobin = [ 0; rho[ii][sp[step:step:nrho]] ]
        jj = findall(.!ismissing.(rhobin) .& isfinite.(rhobin))
        D2inbin = binavg(D2[ii], rho[ii], rhobin[jj])
        rhoinbin = binavg(rho[ii], rho[ii], rhobin[jj])
        nbin_out = min(nbin_out_max, length(rhobin))
        return nbin_out, rhobin[1:nbin_out], D2inbin[1:nbin_out], rhoinbin[1:nbin_out]
    else
        return 1, [missing], [missing], [missing]
    end
end

"""
structure function D2, rho, A, epsilon at each level from w stare
D2bin, rhobin, A, noise = D2_rho_stare( w, pitch, roll, Ur, Vr; out=17 )
"""
function D2_rho_stare( w, pitch, roll, Ur, Vr; nbin_out_max=17 )

    nbin_out = nbin_out_max
    
    (nt, nz) = size(w)
    A      = Vector{Union{Missing,Float64}}(missing, nz)
    noise  = Vector{Union{Missing,Float64}}(missing, nz)
    rhobin = Matrix{Union{Missing,Float64}}(missing, nbin_out, nz)
    D2bin  = Matrix{Union{Missing,Float64}}(missing, nbin_out, nz)
    for izo in 1:nz # loop vertically
        #=
        ci1,ci2, li1,li2, it1,iz1,it2,iz2 = lidarindices(nt, nz, izo) # might do outside the loop
        zm, dr2, dz2, D2 = displacements( ci1,ci2, Ur*timestep,Vr*timestep,
                                          pitch,roll, w; timestep=timestep )
        rho = rhopair.(dr2, dz2) # approx r^2/3
        # bin average str fcn D2 in equally-populated bins of rho
        @show size(rho), size(D2)
        rhobin_, D2inbin_, rhoinbin_ = equal_bin(rho, D2)
        rhobin[:,izo] .= rhoinbin_[1:nbin_out]
        D2bin[ :,izo] .= D2inbin_[ 1:nbin_out]
        # regress to get A
        ii = .!ismissing.(rhobin[:,izo]) .& .!ismissing.(D2bin[:,izo])
        if sum(ii) > 2
            A[izo] = anom(rhobin[:,izo][ii]) \ anom(D2bin[:,izo][ii])
            noise[izo] = mean(D2bin[:,izo][ii]) - A[izo] * mean(rhobin[:,izo][ii]) # noise
        end
        =#
        ci1,ci2, li1,li2, it1,iz1,it2,iz2 = lidarindices(nt, nz, izo) # might do outside the loop
        zm, dr2, dz2, D2 = displacements( ci1,ci2, Ur*timestep,Vr*timestep,
                                          pitch,roll, w; timestep=timestep )
        rho = rhopair.(dr2, dz2) # approx r^2/3
        # bin average str fcn D2 in equally-populated bins of rho
        nbin_actual, rhobin_, D2inbin_, rhoinbin_ = equal_bin(rho, D2; nbin_out_max=nbin_out_max)
        rhobin[1:nbin_actual,izo] .= rhoinbin_
        D2bin[ 1:nbin_actual,izo] .= D2inbin_
        # regress to get A
        ii = .!ismissing.(rhobin[1:nbin_actual,izo]) .& .!ismissing.(D2bin[1:nbin_actual,izo])
        if sum(ii) > 2
            A[izo] = anom(rhobin[1:nbin_actual,izo][ii]) \ anom(D2bin[1:nbin_actual,izo][ii])
            noise[izo] = mean(D2bin[1:nbin_actual,izo][ii]) - A[izo] * mean(rhobin[1:nbin_actual,izo][ii]) # noise
        end
    end
    return D2bin, rhobin, A, noise
end

# functions for subsetting and finding the offset with max covariance
# newer: 2025-02

"""
return indices jl, js to subset windows dtl[jl], dts[js] st. 
limdtl[1]+offset <= dtl[jl] <= limdtl[2]+offset
limdtl[1]        <= dts[js] <= limdtl[2]
offset shifts the window in the long step coordinates
"""
# use with code for chunks found elsewhere
function offset_subset(dtl, dts, limdtl, offset=eltype(limdtl)(0))
    # index the long data set (gappy Halo) with absolute time deltas
    jl = findall(limdtl[1] .<= dtl-offset .<= limdtl[2])
    # comb the time indices out of VN dts
    js = findindices(dtl[jl] .- offset, dts) # findindices( needles, haystack )
    return jl, js
end

function offset_cov(dtl, dts, limdtl, offset, yl, ys)
    jl, js = offset_subset(dtl, dts, limdtl, offset)
    # try
    #     ii = isfinite.(yl[jl]) .&& isfinite.(ys[js]) # skip NaNs # sometimes breaks: arrays not broadcast to consistent size
    #     yl[jl][ii]
    # catch
    #     print("limdtl=$(limdtl) offset=$(offset) jl($(length(jl))), js($(length(js)))")
    # end
    # return cov = mean( skipmissing(yl[jl][ii] .* ys[js][ii]) ) # skip missings
    nn = @. good(yl[jl]) & good(ys[js])
    a_cov = cov(yl[jl][nn], ys[js][nn])
    # a_cov = mean( skipmissing(anom(yl[jl]) .* anom(ys[js])) ) # skip missings
    return a_cov
end

"find optimal offset timedelta (seconds) that syncs yl, ys"
function offset_range_covs(dtl, dts, limdtl, rangeoffset, yl, ys)
    covs = [ offset_cov(dtl, dts, limdtl, offset, yl, ys) 
                for offset in rangeoffset ]
end

"""
return the time offset that syncs dtl (Halo) and dts (VN) for the 
window limdt=[chunkdtstart, chunkdtend]
"""
function sync_offset(dtl, dts, yl, ys, limdt, rangeoffset=Second(0):Second(1):Second(200))
    rangecovs = offset_range_covs(dtl, dts, limdt, rangeoffset, yl, ys)
    maxcov, fm = findmax(rangecovs)
    bestoffset = rangeoffset[fm]
    return bestoffset, maxcov, std(rangecovs)
end

# load input data for doing the offset calculations
"load VN and 1 day of lidar timeangles"
function load_vn_lidar_data(thisdt, Vn=read_vecnav_dict())
    dts = Vn[:vndt] # short timesteps
    ys = Vn[:Roll]  # short-step data

    # read the lidar time axis for this day
    dtstamp = Dates.format(thisdt, dateformat"yyyymmdd")
    datapath = joinpath.(pwd(),"data",dtstamp)
    files = filter(startswith("Stare"), readdir(datapath))
    fullfiles = joinpath.(datapath, files)
    # also read the first hour of the next day
    nextdt = thisdt + Day(1)
    nextdatapath = joinpath.(pwd(), "data", Dates.format(nextdt, "yyyymmdd"))
    if isdir(nextdatapath)
        hour00 = readdir( nextdatapath ) |> filter(startswith("Stare")) |> filter(endswith("_00.hpl"))
        full25files = [fullfiles ; joinpath(nextdatapath, hour00[1])]
    else
        full25files = fullfiles
    end
    # read all times in those files
    ta, _, _ = read_lidar.read_streamlinexr_beam_timeangles(full25files)
    tatime = ta[:time] # lidar time axis, hours
    i20 = findfirst(tatime .> 20.0)
    wrap = (i20-1) .+ findall( tatime[i20:end] .< 5 )
    tatime[wrap] .+= 24.0 # increment wrapped times from next day by 24 h
    dtl = @. thisdt + Millisecond(round(Int64, tatime * 3_600_000)) # long timesteps
    yl = ta[:pitch] # long-step data
    return dtl, yl, dts, ys
end

# read VectorNav
Vn = read_vecnav_dict()

# procedural test of windowing and timing offsets 
rangeoffset = -Second(0):Second(1):Second(200)

thisdt = DateTime(2024,6,8) # select day
dtl, yl, dts, ys = load_vn_lidar_data(thisdt, Vn)

limdt = thisdt + Hour(1) + Minute(1) .+ Minute.([0, 3])

# test computing one covariance
#=
offset = Second(126)

jl = findall(limdt[1] .<= dtl-offset .<= limdt[2])
# comb the time indices out of VN dts
js = findindices(dtl[jl] - offset, dts) # findindices( needles, haystack )
nn = @. good(yl[jl]) & good(ys[js])
a_cov = cov(yl[jl][nn], ys[js][nn])

# test computing a range of covariances
rangecovs = Vector{Float64}(undef, length(rangeoffset))
for (i, offset) in enumerate(rangeoffset)
    jl = findall(limdt[1] .<= dtl-offset .<= limdt[2])
    # comb the time indices out of VN dts
    js = findindices(dtl[jl] - offset, dts) # findindices( needles, haystack )
    nn = @. good(yl[jl]) & good(ys[js])
    rangecovs[i] = cov(yl[jl][nn], ys[js][nn])
end
maxcov, imax = findmax(rangecovs)

clf()
subplot(3,1,1)
plot(Dates.value.(rangeoffset), rangecovs)
plot(Dates.value.(rangeoffset[imax]), maxcov, marker="o")
title("max=$(maxcov), std=$(std(rangecovs))")
gcf()
=#

#   Positive offsets make the l window select from forward in the original l
#   timeseries, and shift the data in this window backward to compare with an
#   earlier (0-offset) time in the s series. Thus l signals are advanced to
#   earlier times to align with earlier (0 offset). signals in the s series.
# 
#   This disagram depicts the alignment of the windows of the Halo and VectorNav
#   data windows.
# 
#                    ---> time
#   Halo    l   ---|----signals |ere----        Halo time series with |requested window|
#               --->>>>[signals here]---   >>>> offset
#               ---[signals here]<<<<---   <<<< aligns later l data with s data
#               ---[vecnavs gnal]-------        nearest points to nonmissing data in window
#   VecNav  s   ---|vecnavsignal|-------

# get all the files, and all the unique hours of the files
allstarefiles = vcat( [ joinpath.("data",F, 
    filter( startswith(r"Stare_"), readdir(joinpath("data",F)) ) ) 
  for F in filter( startswith(r"20240"), readdir("data") ) ]... )

REm = match.(r"Stare_116_(\d{8}_\d{2}).hpl", allstarefiles)
dth = [ DateTime(r[1], dateformat"yyyymmdd_HH") for r in REm ]
unique(floor.(dth, Hour)) # all 991 are already unique

# also linked in ./data/all

# get all times and a main list of chunks
# regardless of Stare file boundaries

# indices of gaps to find stare chunks
function stare.all_gaps(dt::Vector{DateTime})
    ien = findall( diff(dt) .> Second(36) )
    ist = ien .+ 1
    return ien, ist
end

if isfile("lidar_dt.jld2") 
    LidarDt = load("lidar_dt.jld2")
    # Dict: DateTime of each beam and indices of the start and end of staring chunks.
else
    # load and catenate time from the many Stare hpl files
    # construct and save a time vector of datetimes for all lidar beams
    # moved to save_lidar_dt.jl
    starefiles = filter(startswith("Stare_116_"), readdir(joinpath(lidarstemdir, "all")))
    ff = joinpath.(lidarstemdir, "all", starefiles)
    ta, _,_ = read_lidar.read_streamlinexr_beam_timeangles(ff) # slow
    
    day0 = floor(ta[:start_time][1], Day)
    dtime = day0 .+ Millisecond.(round.(Int64, 3600_000 .* ta[:time]))
    for j in findall(diff(ta[:time]) .< -1.0)
        dtime[(j+1):end] .+= Day(1) # increment days
    end

    ien, ist = all_gaps(dtime)
    ist = [1; ist]
    ien = [ien; length(dtime)]
    # [ist ien]
    
    @save "lidar_dt.jld2" dtime ist ien
    LidarDt = load("lidar_dt.jld2")
end

# part out the data among the individual files
lidarstemdir = "./data"
starefiles = filter(startswith("Stare_116_"), readdir(joinpath(lidarstemdir, "all")))
ff = joinpath.(lidarstemdir, "all", starefiles)

if isfile("file_beam_inds.jld2")
    FileInds = load("file_beam_inds.jld2")
else # slow reload
    nfiles = length(ff)
    nbeams = zeros(Int32, nfiles)
    ngates = -1
    nheaderlines = 17
    # read number of lines for each file
    for (i,file) in enumerate(ff)
        global h, ngates
        h = read_lidar.read_streamlinexr_head(file)
        nlines = h[:nlines]
        ngates = h[:ngates]
        # beams could be rays or times
        nbeams[i] = round(Int, (nlines - nheaderlines) / (1+ngates)) # number of beams for each file
    end

    # indices for files slices
    bigind_file_end   = cumsum(nbeams)
    bigind_file_start = [0; bigind_file_end[1:end-1]] .+ 1

    # times, ibeam1, nbeams, hdr = read_lidar.read_streamlinexr_beam_times(ff) # slow bc. it rereads times
    # all( bigind_file_start .== ibeam1 ) # should be true
    [bigind_file_start nbeams bigind_file_end]
    @save "file_beam_inds.jld2" bigind_file_start nbeams bigind_file_end
    FileInds = load("file_beam_inds.jld2")
end

# load indices
LidarDt = load("lidar_dt.jld2")
length(LidarDt["dtime"]) == FileInds["bigind_file_end"][end] # should be true

# define periodic data types

# Wrap indices for periodic behavior
_wrap(x::Integer, n::Int) = mod1(x, n)
_wrap(x::AbstractRange{<:Integer}, n::Int) = mod1.(collect(x), n)
_wrap(x::AbstractArray{<:Integer}, n::Int) = mod1.(x, n)
_wrap(::Colon, ::Int) = (:)

# PeriodicVector periodically indexes
struct PeriodicVector{T}
    data::Vector{T}
end

# convenience constructor
PeriodicVector{T}(x, n::Integer) where {T} =
    PeriodicVector{T}(fill(x, n))

Base.ndims(p::PeriodicVector) = ndims(p.data) # literally 1
Base.size(p::PeriodicVector) = size(p.data)
Base.length(p::PeriodicVector) = length(p.data)
Base.eltype(::Type{PeriodicVector{T}}) where {T} = T

Base.getindex(p::PeriodicVector, i) =
    getindex(p.data, _wrap(i, length(p.data)))

Base.setindex!(p::PeriodicVector, v, i) =
    setindex!(p.data, v, _wrap(i, length(p.data)))
Base.setindex!(p::PeriodicVector, v::AbstractVector, i::Union{AbstractVector, AbstractRange}) =
    (p.data[_wrap(i, length(p.data))] .= v)
    
Base.iterate(P::PeriodicVector, state...) = iterate(P.data, state...)
Base.show(io::IO, P::PeriodicVector) = show(io, P.data)

# PeriodicMatrix periodically indexes the first dimension
struct PeriodicMatrix{T}
    data::Matrix{T}
end

PeriodicMatrix{T}(x, nrows::Integer, ncols::Integer) where {T} =
    PeriodicMatrix{T}(fill(x, nrows, ncols))

Base.ndims(p::PeriodicMatrix) = ndims(p.data) # literally 2
Base.size(P::PeriodicMatrix) = size(P.data)
Base.eltype(::Type{PeriodicMatrix{T}}) where {T} = T

Base.getindex(P::PeriodicMatrix, i, j...) =
    getindex(P.data, _wrap(i, size(P.data,1)), j...)

Base.setindex!(P::PeriodicMatrix, v, i, j...) =
    setindex!(P.data, v, _wrap(i, size(P.data,1)), j...)
Base.setindex!(P::PeriodicMatrix, v::AbstractArray, i::Union{AbstractVector, AbstractRange}, j...) =
    (P.data[_wrap(i, size(P.data,1)), j...] .= v)

    # add iteration, show
Base.iterate(P::PeriodicMatrix, state...) = iterate(P.data, state...)
Base.show(io::IO, P::PeriodicMatrix) = show(io, P.data)


"initialize a beams Dict with exactly nbeams beams with periodic 1st (time) index"
function init_periodic_beams(nbeams, ngates)
    # helper functions to initialize periodic data types
    PeriodicVector_missing(n) = PeriodicVector(
        Vector{Union{Float32,Missing}}(fill(missing, n)) )
    PeriodicMatrix_missing(nrows, ncols) = PeriodicMatrix(
        Matrix{Union{Float32,Missing}}(fill(missing, nrows, ncols)) )
    Vector_missing(n) = Vector{Union{Float32,Missing}}(fill(missing, n))

    beams = Dict(
        :time       => PeriodicVector_missing(nbeams), # decimal hours
        :azimuth    => PeriodicVector_missing(nbeams), # degrees
        :elevangle  => PeriodicVector_missing(nbeams),
        :pitch      => PeriodicVector_missing(nbeams),
        :roll       => PeriodicVector_missing(nbeams),
        :height     =>         Vector_missing(ngates), # center of gate
        :dopplervel => PeriodicMatrix_missing(nbeams,ngates), # m/s
        :intensity  => PeriodicMatrix_missing(nbeams,ngates), # SNR + 1
        :beta       => PeriodicMatrix_missing(nbeams,ngates) # m-1 sr-1  backscatter?
        )
end

# periodic buffer loop for data
nx = 4000 
nz = 80
x = zeros(nx,nz) # Doppler vel, backscatter, etc.
P = PeriodicMatrix( x )
dtx = PeriodicVector( fill(DateTime(0), nx) )

beams = init_periodic_beams(nx, nz)
beams[:time] isa Type # should be an instance, not a type

"""
modifying read_streamlinexr_stare!(file_path, header, beams, bb)
Read data and fill in the beams for a single file.
"""
function read_streamlinexr_stare!(file_path, h, beams, bb, nheaderlines=17; startat=1, endat=0)
    # beams is a Dict of PeriodicVector and PeriodicMatrix
    # bb is the big_index range, which will be interpreted periodically by arrays in beams

    # use header information in h
    nz = size(beams[:height][:],1)

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

    # parse the variables into the dict beams by beam
    # beams[:time     ][bb] .= beam_timeangles[:,1] # decimal hours
    # beams[:azimuth  ][bb] .= beam_timeangles[:,2] # degrees
    # beams[:elevangle][bb] .= beam_timeangles[:,3] # degrees
    # beams[:pitch    ][bb] .= beam_timeangles[:,4]
    # beams[:roll     ][bb] .= beam_timeangles[:,5]
    # # by gate
    # beams[:height   ][1:nz] .= (beam_velrad[1,1:nz,1].+0.5) .* h[:gatelength] # center of gate, assumes same for all beams

    # # dependent variables (beam, gate)
    # beams[:dopplervel][bb,1:nz] .= beam_velrad[:,1:nz,2] # m/s
    # beams[:intensity ][bb,1:nz] .= beam_velrad[:,1:nz,3] # SNR + 1
    # beams[:beta      ][bb,1:nz] .= beam_velrad[:,1:nz,4] # m-1 sr-1  backscatter

    setindex!(beams[:time],      beam_timeangles[:,1], bb) # decimal hours
    setindex!(beams[:azimuth],   beam_timeangles[:,2], bb) # degrees
    setindex!(beams[:elevangle], beam_timeangles[:,3], bb) # degrees
    setindex!(beams[:pitch],     beam_timeangles[:,4], bb) 
    setindex!(beams[:roll],      beam_timeangles[:,5], bb)
    # by gate
    beams[:height][1:nz] .= (beam_velrad[1,1:nz,1].+0.5) .* h[:gatelength] # center of gate
    setindex!(beams[:dopplervel], beam_velrad[:,1:nz,2], bb, 1:nz) # m/s
    setindex!(beams[:intensity],  beam_velrad[:,1:nz,3], bb, 1:nz) # SNR + 1
    setindex!(beams[:beta],       beam_velrad[:,1:nz,4], bb, 1:nz) # m-1 sr-1  backscatter
end

nannoise(noisemask, vel) = ismissing(noisemask) ? missing : (noisemask ? vel : missing)

"plot backscatter and velocity"
function pcolor_lidar_stare(fig, beams, LidarDt, st, en, noisethr=1.03)
    dt = LidarDt["dtime"][st]
    dstr = Dates.format(dt, dateformat"yyyymmdd_HHMM")
    
    # get data
    height = beams[:height]
    # subset the variables for the present chunk
    tmp = beams[:time][st:en] # decimal hours
    # re-center around first time
    time = tmp[1] .+ (mod.(tmp .- tmp[1] .+ 0.5, 1.0) .- 0.5) # hours
    kys = Symbol.(split("beta dopplervel intensity"))
    (beta, dopplervel, intensity) = (beams[k][st:en,:] for k in kys)
    noisemask = intensity .>= noisethr
    
    # remove the vertical-mean Doppler velocity
    # mdv = mapslices(x->mean(skipmissing(x)), dopplervel; dims=2)[:]
    # vel = dopplervel .- mdv # subtract mean heave vel
    # noisemask = intensity.>=noisethr
    f(x) = !ismissing(x) && isfinite(x) && x>noisethr
    vel = remove_mdv(f, intensity, dopplervel)

    clf()
    # plot_stare(time, height, beta, dopplervel, intensity)
    subplot(2,1,1)
    pcolormesh(time*60, height[4:end]/1e3, pd(beta[:,4:end]),
        cmap=ColorMap("RdYlBu_r"), vmin=0, vmax=0.2*maximum(beta[:,4:end]))
    xtl = round(time[1]*60) : 1 : round(time[end]*60)
    xticks(xtl)
    gca().set_xticklabels(xtl.%60)
    ylim([0, 2])
    ylabel("height (km)")
    xlabel("time (minute)")
    title("$(Dates.format(dt, dateformat"yyyy-mm-dd HH:MM"))\nbackscatter")
    colorbar()

        # true for valid, false for noise
    subplot(2,1,2)
    pcolormesh(time*60, height[4:end]/1e3, pd(nannoise.(noisemask, vel)[:,4:end]),
        cmap=ColorMap("RdYlBu_r"), vmin=-2, vmax=2)
    xtl = round(time[1]*60) : 1 : round(time[end]*60)
    xticks(xtl)
    gca().set_xticklabels(xtl.%60)
    ylim([0, 2])
    ylabel("height (km)")
    xlabel("time (minute)")
    title("Doppler velocity (m/s)")
    colorbar()
    
    tight_layout()
end

 # check indices loop through chunks. advance to next file as data is needed
 # BUT do nothing
#=
let
    bigind_file_ends = FileInds["bigind_file_end"]
    bigind_file_starts  = FileInds["bigind_file_start" ]
    iens             = LidarDt[ "ien"            ]
    ists             = LidarDt[ "ist"            ]

    ifile = 0
    bigind_file_end = 0 # forces initial read in loop
    # bigind_file_end = bigind_file_ends[ifile]
    # bigind_file_st  = bigind_file_starts[ ifile]
    # load the file ff[ifile] data
    for ic in eachindex(iens) # loop over all chunks in the record
        ien = iens[ic]; ist = ists[ic]
        if ien > bigind_file_end
            # increment file
            ifile += 1
            bigind_file_end = bigind_file_ends[ifile]
            bigind_file_start  = bigind_file_starts[ ifile]
            print("\nfile $(ifile) $(bigind_file_start)---$(bigind_file_end) chunk ")
            # load the file ff[ifile] data
        end
        # process chunk
        print("$(ist)-$(ien) ")
    end
end
=#

# load cruiselong data
if !@isdefined UV # | true # do once, otherwise save time
    # load all-2024 relative horizontal winds
    UV = NCDataset(joinpath("data/netcdf", "ekamsat_lidar_uv_20240428-20240604.nc")) # NCDataset

    bigind_file_ends    = FileInds["bigind_file_end"]
    bigind_file_starts  = FileInds["bigind_file_start"]
    iens             = LidarDt[ "ien"            ]
    ists             = LidarDt[ "ist"            ]
end

LidarDt["dtime"]
dtime_st = LidarDt["dtime"][ists]
dtime_en = LidarDt["dtime"][iens]

# chunk indices for which VectorNav data is available
icvn = findfirst( dtime_st .>= Vn[:vndt][1] ):findlast( dtime_en .<= Vn[:vndt][end] )

# loop through all chunks, load next file as data is needed

# periodic buffer for data
nx = 4000 
nz = 80

# Our way of looping does not synchronize chunks to files a priori
# but loads files as needed. 
# queue_indices queues up indices for record sample index i:
"ic, ifile, ist, ien, bigind_file_start, bigind_file_end = queue_indices(i)"
function queue_indices(i, 
    ists=ists, iens=iens, 
    bigind_file_starts=bigind_file_starts, bigind_file_ends=bigind_file_ends)
    ic = findfirst(iens .>= i)
    ifile = findfirst(bigind_file_ends .>= i)
    bigind_file_start = bigind_file_starts[ifile]
    bigind_file_end   = bigind_file_ends[  ifile]
    ic, ifile, ists[ic], iens[ic], bigind_file_start, bigind_file_end
end

ic, ifile, ist, ien, bigind_file_start, bigind_file_end = queue_indices(2919682)
# call queue_indices here to start in mid stream.

let # local scope
    # x = zeros(nx,nz) # Doppler vel, backscatter, etc.
    # P = PeriodicMatrix( x )
    # dtx = PeriodicVector( fill(DateTime(0), nx) )
    h = read_lidar.read_streamlinexr_head(ff[1]) # initialize globals
    beams = init_periodic_beams(nx, nz)

    # TKE dissipation chunk output data array
    epsi_tmp = zeros(Float64, length(iens), nz) .- 5 
    # -5 is uncomputed missing value.
    # Record posibilities for why dissipation could be missing as sentinel
    # negative values, since physical TKE dissipation is always positive.
    indmiss(x) = ismissing(x) ? -3 : x # individually missing calculated values = -3

    fig = gcf()

    # start from beginning of cruise
    ifile = 1
    bigind_file_end = bigind_file_ends[ifile] # forces initial read in loop
    bigind_file_start = bigind_file_starts[ifile]
    
    # queue up a particular time
    # ic, ifile, ist, ien, bigind_file_start, bigind_file_end = queue_indices(2919682)
    # ic, ifile, ist, ien, bigind_file_start, bigind_file_end = queue_indices(3094111) # last file

    #=
    # test nasssty litte casessses
    # special treatment for starting in the middle
    ifile = 187
    bigind_file_end = bigind_file_ends[ifile]
    bigind_file_start = bigind_file_starts[ifile]
    =#
    # initial read
    tmp = read_lidar.read_streamlinexr_head(ff[ifile])
    [ h[k] = tmp[k] for k in keys(h) ] # update h wo rebinding it
    bb = bigind_file_start:bigind_file_end
    read_streamlinexr_stare!(ff[ifile], h, beams, bb)
    for ic in eachindex(iens) # ic:icvn[end] # 1118:5793 # loop over all chunks in the record
        ien = iens[ic]; ist = ists[ic]
        if ien > bigind_file_end # need to load more data; if doesn't create local scope
            while ien > bigind_file_end # need to load more data
                ifile += 1
                bigind_file_end = bigind_file_ends[ifile]
            end
            bigind_file_start  = bigind_file_starts[ifile]
            print("\nfile $(ifile) $(bigind_file_start)---$(bigind_file_end) chunk ")
            # load data
            tmp = read_lidar.read_streamlinexr_head(ff[ifile])
            [ h[k] = tmp[k] for k in keys(h) ] # update h wo rebinding it
            bb = bigind_file_start:bigind_file_end
            read_streamlinexr_stare!(ff[ifile], h, beams, bb)
            # beams[:key][bb] is the data loaded for this file
            # beams[:key][ist:ien] is the data for a chunk
        end
        # process one chunk
        print("$(ist)-$(ien) ")

        # compute dissipation for the chunk
        # try # read a chunk, collocate wind and VN data
            # read_stare_chunk organizes and aligns lidar, motion, and wind data
            # should also work for periodic arrays in beams
            dt = Date(h[:start_time])
            # print(beams[:time][ist], ", ",beams[:time][ien], ", ") # OK
            # print(UV[:ur][1,1]) # probably OK
            dopplervel, pitch, roll, vn0, vn1, vn2, Ur, Vr, mdv = read_stare_chunk( dt, beams, Vn, UV, ist, ien )
            if any(isfinite.(Ur)) && any(isfinite.(Vr)) # there is wind data
                w = dopplervel .- mdv

                len = minimum([size(w,1), size(pitch,1), size(Ur,1)])  # trim to shortest length
                D2bin, rhobin, A, noise = D2_rho_stare( w[1:len,:], pitch[1:len]*pi/180, roll[1:len]*pi/180, Ur[1:len,:], Vr[1:len,:] )
                epsi_tmp[ic,:] .= @. indmiss( epsilon(max(0,A)) )
                # sets individually missing dissipation values to -3
            else
                epsi_tmp[ic,:] .= -4 # code for missing wind
            end
        # catch
        #         epsi_tmp[ic,:] .= -5 # code for missing data, probably VN missing
        # end

        # make and save figure for the chunk, stamped with the chunk start time
        #=
        clf()
        pcolor_lidar_stare(fig, beams, LidarDt, ist, ien)
        dstr = Dates.format(LidarDt["dtime"][ist], dateformat"yyyymmdd_HHMM")
        savefig(joinpath("data","plot","Chunk1_$(dstr).png"))
        =#
    end

    epsi = epsi_tmp
    height = beams[:height]
    savetime = Dates.format(now(), dateformat"yyyymmdd_HHMM")
    @save joinpath("epsilon_data","epsi_stare_chunks_r$(savetime).jld") epsi dtime_st dtime_en height
    print(savetime)

    clf()
    pcolormesh(dtime_st[1:16], height, pd(epsi[1:15,1:79]))
    gcf()

end

# test
#=
test = load("epsilon_data/epsi_stare_chunks_r$(savetime).jld")

norm = PyPlot.matplotlib.colors.LogNorm(vmin=1e-5, vmax=1e-2)
clf()
pcolormesh(test["dtime_en"], test["height"], pd(max.(1e-6, test["epsi"])), norm=norm, cmap=ColorMap("RdYlBu_r"))
gcf()
=#
