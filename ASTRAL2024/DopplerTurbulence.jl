module DopplerTurbulence

using Revise
using Pkg; Pkg.activate(".")

using Dates
using Statistics
using Rotations
using Interpolations
using DSP
using FFTW

# Include and import uncertainty propagation functions
include("displacement_uncertainty.jl")
using .DisplacementUncertainty

# Export DisplacementUncertainty functions for use in notebooks
export displacement_variance_pitch_roll, propagate_rho_uncertainty,
       combine_variances, estimate_pitch_roll_uncertainty

# Export utility functions for use in notebooks
export pd, m2n, n2m, missmean, anom, binavg, hp, findindices, indavg,
       trigs, wtrue_trigs, wtrue, uniquepairs, allcross, rng, lidarindices,
       rangegate, displacements, rhopair
export epsilon, fit_valid_xy
export epsilon_ci95_from_a_ci, fit_stats_onepass, trim_structure_inputs, equal_bin
export tls_slope_intercept, D2_rho_stare, D2_rho_stare_qc
export small_displacements, small_bin, D2_rho_small_disp
export transfer_fcn, get_noise, get_slope, get_epsilon

# utility functions
pd = permutedims
m2n(x) = ismissing(x) ? NaN : x

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
        bi = searchsortedlast(b, x)
        bi = clamp(bi, 1, length(b))
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

"average xl within windows to right of points of the index ind of xl"
function indavg(xl, ind; full=20)
    xm = zeros(Float64, size(ind))
    for (i,idx) in enumerate(ind)
        ii = max(1,idx) : min(length(xl),idx+full)
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

# displacements with no adjustment for tilting into the horizontal wind 
# U, V vary slowly; pitch,roll,w vary fast
# there are nt*(nt-1)/2 ~ O(nt^2) outputs, so correct stuff first


# functions for structure functions

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

# displacments and structure functions 

rangegate = 24.0 # for ASTRAL 2024 Halo Photonics StreamLineXR

"""
zm, dr2, dz2, D2 = displacements(ci1,ci2,it1,iz1,it2,iz2,Udt,Vdt,pitch,roll,w; rangegate=rangegate)
Displacements of sample pairs for one (vertical) subvolume using precomputed pair-index vectors.
"""
function displacements(ci1, ci2, it1, iz1, it2, iz2, Udt, Vdt, pitch, roll, w; rangegate=rangegate)
    rng(iz) = rangegate * (iz - 1 + 0.5)

    # horiz translation of the sample volumes by mean wind
    Udtbar = @. (Udt[iz2] + Udt[iz1]) / 2
    Vdtbar = @. (Vdt[iz2] + Vdt[iz1]) / 2
    X = @. Udtbar * (it2 - it1)
    Y = @. Vdtbar * (it2 - it1)

    # vertical middle of pair
    zm = @. (rng(iz2) * cos(pitch[it2]) * cos(roll[it2]) + rng(iz1) * cos(pitch[it1]) * cos(roll[it1])) / 2

    # displacement between pair of points
    dz = @. rng(iz2) * cos(pitch[it2]) * cos(roll[it2]) - rng(iz1) * cos(pitch[it1]) * cos(roll[it1])
    dx = @. X + rng(iz2) * -sin(pitch[it2]) - rng(iz1) * -sin(pitch[it1])
    dy = @. Y + rng(iz2) * cos(pitch[it2]) * sin(roll[it2]) - rng(iz1) * cos(pitch[it1]) * sin(roll[it1])

    dz2 = dz .* dz
    dr2 = @. dz2 + dx * dx + dy * dy
    D2 = @. (w[ci2] - w[ci1])^2
    return zm, dr2, dz2, D2
end

"""
zm, dr2, dz2, D2, var_dr2, var_dz2 = displacements(ci1,ci2,it1,iz1,it2,iz2,Udt,Vdt,pitch,roll,w; σ_pitch=σ_pitch, σ_roll=σ_roll, rangegate=rangegate)
Computes structure function D2 for
Displacements of sample pairs for one (vertical) subvolume using precomputed pair-index vectors,
including variance estimates.
"""
function displacements(ci1, ci2, it1, iz1, it2, iz2,
            Udt, Vdt, pitch, roll, w,
            σ_pitch::Real, σ_roll::Real; rangegate=rangegate)
    rng(iz) = rangegate * (iz - 1 + 0.5)
    Udtbar = @. (Udt[iz2] + Udt[iz1]) / 2
    Vdtbar = @. (Vdt[iz2] + Vdt[iz1]) / 2
    X = @. Udtbar * (it2 - it1)
    Y = @. Vdtbar * (it2 - it1)

    zm = @. (rng(iz2) * cos(pitch[it2])*cos(roll[it2]) + rng(iz1) * cos(pitch[it1])*cos(roll[it1])) / 2

    dz = @. rng(iz2) * cos(pitch[it2])*cos(roll[it2]) - rng(iz1) * cos(pitch[it1])*cos(roll[it1])
    dx = @. X + rng(iz2) *-sin(pitch[it2]) - rng(iz1) *-sin(pitch[it1])
    dy = @. Y + rng(iz2) * cos(pitch[it2])*sin(roll[it2]) - rng(iz1) * cos(pitch[it1])*sin(roll[it1])

    dz2 = dz .* dz
    dr2 = @. dz2 + dx*dx + dy*dy
    D2 = @. (w[ci2] - w[ci1])^2

    # NEW: Compute uncertainty contributions from pitch/roll
    var_dr2 = similar(dr2)
    var_dz2 = similar(dz2)

    for i in eachindex(dr2)
        var_dr2[i], var_dz2[i] = displacement_variance_pitch_roll(
            rng(iz1[i]), rng(iz2[i]),
            pitch[it1[i]], pitch[it2[i]],
            roll[it1[i]], roll[it2[i]],
            σ_pitch, σ_roll
        )
    end

    return zm, dr2, dz2, D2, var_dr2, var_dz2
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

"return finite non-missing x,y used in linear fit"
function fit_valid_xy(x, y)
    ii = .!ismissing.(x) .& .!ismissing.(y) .& isfinite.(x) .& isfinite.(y)
    return Float64.(x[ii]), Float64.(y[ii])
end

"propagate A 95% CI to epsilon 95% CI via epsilon(A) transform"
function epsilon_ci95_from_a_ci(a, alo, ahi)
    if ismissing(a) || ismissing(alo) || ismissing(ahi)
        return missing, missing
    end
    if !(isfinite(a) && isfinite(alo) && isfinite(ahi))
        return missing, missing
    end
    if a <= 0
        return missing, missing
    end

    lo = min(alo, ahi)
    hi = max(alo, ahi)
    lo_pos = max(lo, 0.0)
    hi_pos = max(hi, 0.0)
    if hi_pos <= 0
        return missing, missing
    end

    return epsilon(lo_pos), epsilon(hi_pos)
end

"single-pass slope and quality statistics from binned x,y"
function fit_stats_onepass(x, y)
    xv, yv = fit_valid_xy(x, y)
    nbins = length(xv)

    out_missing = (
        A=missing,
        noise=missing,
        R2=missing,
        RMSE=missing,
        se_A=missing,
        A_ci_lo=missing,
        A_ci_hi=missing,
        epsi_ci_lo=missing,
        epsi_ci_hi=missing,
        nbins=nbins,
    )

    if nbins < 3
        return out_missing
    end

    xbar = mean(xv)
    ybar = mean(yv)
    xa = xv .- xbar
    ya = yv .- ybar

    sxx = sum(abs2, xa)
    sxy = sum(xa .* ya)
    sst = sum(abs2, ya)
    if !isfinite(sxx) || !isfinite(sxy) || !isfinite(sst) || sxx <= 0 || sst <= 0
        return out_missing
    end

    A = sxy / sxx
    if !isfinite(A)
        return out_missing
    end
    noise = ybar - A * xbar

    resid = ya .- A .* xa
    sse = sum(abs2, resid)
    dof = nbins - 2
    if !isfinite(sse) || sse < 0 || dof <= 0
        return (; out_missing..., A=A, noise=noise)
    end

    mse = sse / dof
    if !isfinite(mse) || mse < 0
        return (; out_missing..., A=A, noise=noise)
    end

    R2 = 1 - sse / sst
    RMSE = sqrt(mse)
    se_A = sqrt(mse / sxx)
    if !isfinite(RMSE) || !isfinite(se_A)
        return (; out_missing..., A=A, noise=noise, R2=R2)
    end

    tcrit = quantile(TDist(dof), 0.975)
    delta = tcrit * se_A
    A_ci_lo = A - delta
    A_ci_hi = A + delta
    epsi_ci_lo, epsi_ci_hi = epsilon_ci95_from_a_ci(A, A_ci_lo, A_ci_hi)

    return (
        A=A,
        noise=noise,
        R2=R2,
        RMSE=RMSE,
        se_A=se_A,
        A_ci_lo=A_ci_lo,
        A_ci_hi=A_ci_hi,
        epsi_ci_lo=epsi_ci_lo,
        epsi_ci_hi=epsi_ci_hi,
        nbins=nbins,
    )
end

"trim chunk inputs to the common time length required by structure-function code"
function trim_structure_inputs(w, pitch, roll, Ur, Vr)
    len = minimum((
        size(w, 1),
        length(pitch),
        length(roll),
        size(Ur, 1),
        size(Vr, 1),
    ))
    if len < 60
        return nothing
    end
    return (
        w = w[1:len, :],
        pitch = pitch[1:len],
        roll = roll[1:len],
        Ur = Ur[1:len, :],
        Vr = Vr[1:len, :],
    )
end

"bin average D2 in nbin_out_max equally-populated bins by rho"
function equal_bin(rho, D2; nbin=200, nbin_out_max=17)
    ii = findall(.!ismissing.(rho) .& isfinite.(rho) .& .!ismissing.(D2) .& isfinite.(D2))
    nrho = length(ii) # need to keep track of this total number of realizations
    if nrho >= 8
        sp = sortperm(rho[ii])
        step = max(1, round(Int32, nrho / nbin))
        rhobin = [0; rho[ii][sp[step:step:nrho]]]
        jj = findall(.!ismissing.(rhobin) .& isfinite.(rhobin))
        D2inbin,  D2varbin  = binavgvar(D2[ii], rho[ii], rhobin[jj])
        rhoinbin, rhovarbin = binavgvar(rho[ii], rho[ii], rhobin[jj])
        nbin_out = min(nbin_out_max, length(rhobin))
        return nbin_out, rhobin[1:nbin_out], D2inbin[1:nbin_out], rhoinbin[1:nbin_out], D2varbin[1:nbin_out], rhovarbin[1:nbin_out]
    else
        return 1, [missing], [missing], [missing], [missing], [missing]
    end
end

"equal_bin(rho, D2, rho_var) equal_bin with rho_var: adds propagated instrument uncertainty to rhovarbin"
#  from nbin=200, nbin_out_max=17
function equal_bin(rho, D2, rho_var; nbin=200, nbin_out_max=17)
    rho_var_finite = coalesce.(rho_var, 0.0)  # missing → 0 (no instrument error)
    ii = findall(.!ismissing.(rho) .& isfinite.(rho) .& .!ismissing.(D2) .& isfinite.(D2))
    nrho = length(ii)
    if nrho >= 8
        sp = sortperm(rho[ii])
        step = max(1, round(Int32, nrho / nbin))
        rhobin = [0; rho[ii][sp[step:step:nrho]]]
        jj = findall(.!ismissing.(rhobin) .& isfinite.(rhobin))
        D2inbin,  D2varbin         = binavgvar(D2[ii],              rho[ii], rhobin[jj])
        rhoinbin, rhovarbin_sample = binavgvar(rho[ii],             rho[ii], rhobin[jj])
        rho_var_mean, _            = binavgvar(rho_var_finite[ii],  rho[ii], rhobin[jj])
        nbin_out = min(nbin_out_max, length(rhobin))
        return nbin_out, rhobin[1:nbin_out], D2inbin[1:nbin_out], rhoinbin[1:nbin_out],
               D2varbin[1:nbin_out], (rhovarbin_sample .+ rho_var_mean)[1:nbin_out]
    else
        return 1, [missing], [missing], [missing], [missing], [missing]
    end
end


"""
    binned_tls(r, d, sigma_r, sigma_d)

Calculates the slope and intercept for binned data with errors in both 
variables using Weighted Total Least Squares.
"""
function tls_slope_intercept(r::Vector, d::Vector, var_r::AbstractVector, var_d::AbstractVector; min_r_var=1e-6, min_d_var=1e-6)
    r_bar = mean(r)
    d_bar = mean(d)

    if all(r .== r_bar) || all(d_bar .== d) # singular matrix
        return missing, missing
    end
    
    # A is the r anomaly vector as a N x 1 matrix
    A = reshape(r.-r_bar, :, 1)
    # Variance matrices (Diagonal for independent bins)
    Qaa = Diagonal(max.(var_r, min_r_var))
    Qyy = Diagonal(max.(var_d, min_d_var))
    
    # Cross-covariance is zero for independent r and d errors
    n = length(d)
    Qay = spzeros(n, n)
    
    # Solve for the slope m
    # weighted TLS wtls returns a vector of coefficients; here it's just [slope]
    m = wtls(A, d.-d_bar, Qaa, Qay, Qyy)[1]
    
    # the intercept
    b = d_bar - m * r_bar
    
    return m, b
end
function tls_slope_intercept(r::Vector, d::Vector, sigma_r::Real, sigma_d::Real) 
    n = length(r)
    tls_slope_intercept(r, d, fill(sigma_r, n), fill(sigma_d, n))
end

"""
structure function D2, rho, A, epsilon at each level from w stare
D2bin, rhobin, A, noise = D2_rho_stare(w, pitch, roll, Ur, Vr; out=17)
"""
function D2_rho_stare(w, pitch, roll, Ur, Vr; nbin_out_max=17, method="tls")
    trimmed = trim_structure_inputs(w, pitch, roll, Ur, Vr)
    nbin_out = nbin_out_max
    out_D2 = Matrix{Union{Missing, Float64}}(missing, nbin_out, size(w, 2))
    out_rho = Matrix{Union{Missing, Float64}}(missing, nbin_out, size(w, 2))
    out_A = Vector{Union{Missing, Float64}}(missing, size(w, 2))
    out_noise = Vector{Union{Missing, Float64}}(missing, size(w, 2))

    if isnothing(trimmed)
        return out_D2, out_rho, out_A, out_noise
    end

    w = trimmed.w
    pitch = trimmed.pitch
    roll = trimmed.roll
    Ur = trimmed.Ur
    Vr = trimmed.Vr

    (nt, nz) = size(w)
    A      = Vector{Union{Missing, Float64}}(missing, nz)
    noise  = Vector{Union{Missing, Float64}}(missing, nz)
    nrho   = zeros(Int, nz)
    rhobin = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)
    D2bin  = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)
    rhovarbin = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)
    D2varbin  = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)


    # Item 1 and 2: reuse pair indices and pass precomputed index vectors into displacements.
    level_cache = get_sf_level_cache(nt, nz)
    # Item 3: precompute scaled wind once per call.
    Udt = Ur * timestep
    Vdt = Vr * timestep

    for izo in 1:nz
        lc = level_cache[izo]
        rhoscale_max = max(3500.0, 5.0 * rng(izo))^(2/3)
        zm, dr2, dz2, D2 = displacements(lc.ci1, lc.ci2, lc.it1, lc.iz1, lc.it2, lc.iz2, Udt, Vdt, pitch, roll, w)
        rho = rhopair.(dr2, dz2)

        # total number of valid realizations
        nrho[izo] = sum(.!ismissing.(rho) .& isfinite.(rho) .& .!ismissing.(D2) .& isfinite.(D2))

        nbin_actual, rhobin_, D2inbin_, rhoinbin_, D2varbin_, rhovarbin_ = equal_bin(rho, D2; nbin_out_max=nbin_out_max)
        rhobin[1:nbin_actual, izo] .= rhoinbin_
        D2bin[1:nbin_actual, izo] .= D2inbin_
        rhovarbin[1:nbin_actual, izo] .= rhovarbin_
        D2varbin[1:nbin_actual, izo] .= D2varbin_

        # good values to condition matrices properly
        ii = (    .!ismissing.(rhobin[1:nbin_actual, izo])
               .& .!ismissing.(D2bin[1:nbin_actual, izo])
               .& .!ismissing.(rhovarbin[1:nbin_actual, izo])
               .&  ( rhovarbin[1:nbin_actual, izo] .> 0 )
               .& .!ismissing.(D2varbin[1:nbin_actual, izo])
               .&  ( D2varbin[1:nbin_actual, izo] .> 0 )
               )
        if sum(ii) >= 5 # increased from 2
            if method == "ols"
                # ordinary least squares fit of D2 = A * rho + noise
                A[izo] = anom(rhobin[1:nbin_actual, izo][ii]) \ anom(D2bin[1:nbin_actual, izo][ii])
                noise[izo] = mean(D2bin[1:nbin_actual, izo][ii]) - A[izo] * mean(rhobin[1:nbin_actual, izo][ii])
            elseif method == "tls"
                # total least squares fit of D2 = A * rho + noise, with errors in both variables
                try
                    # use 
                    A_, noise_ = tls_slope_intercept(
                                        rhobin[1:nbin_actual, izo][ii], 
                                        D2bin[1:nbin_actual, izo][ii], 
                                        rhovarbin[1:nbin_actual, izo][ii], 
                                        D2varbin[1:nbin_actual, izo][ii] )
                    A[izo] = A_
                    noise[izo] = noise_
                catch e
                    nrmissing = sum(ismissing.(rhovarbin[1:nbin_actual, izo][ii])) # should be 0
                    ndmissing = sum(ismissing.(D2varbin[1:nbin_actual, izo][ii])) # should be 0
                    nrzero = sum((rhovarbin[1:nbin_actual, izo][ii].==0)) # should be 0
                    ndzero = sum((D2varbin[1:nbin_actual, izo][ii].==0)) # should be 0
                    @warn "TLS failed for level $izo with $nrmissing missing $nrzero zero rhovarbin and $ndmissing missing $ndzero zero D2varbin"
                end
            end
        end
    end
    return D2bin, rhobin, A, noise, nrho
end

"""
D2_rho_stare_qc(w, pitch, roll, Ur, Vr, σ_pitch, σ_roll; nbin_out_max=25)
Like D2_rho_stare, but also returns fit-quality metrics and confidence intervals.
No gating is applied; A and epsilon behavior remains unchanged.
"""
function D2_rho_stare_qc(w, pitch, roll, Ur, Vr, 
                         σ_pitch::Real, σ_roll::Real;
                         nbin_out_max=25, method="tls", Lscale_max=3500.0)

    trimmed = trim_structure_inputs(w, pitch, roll, Ur, Vr)
    nbin_out = nbin_out_max
    nz_out = size(w, 2)
    out_A = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_noise = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_nrho = zeros(Int, nz_out)
    out_R2 = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_RMSE = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_nbins = zeros(Int, nz_out)
    out_se_A = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_A_lo = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_A_hi = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_epsi_lo = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_epsi_hi = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_rho = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz_out)
    out_D2 = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz_out)

    if isnothing(trimmed)
        return out_D2, out_rho, out_A, out_noise, out_R2, out_RMSE, out_nbins, out_se_A, out_A_lo, out_A_hi, out_epsi_lo, out_epsi_hi, out_nrho
    end

    w = trimmed.w
    pitch = trimmed.pitch
    roll = trimmed.roll
    Ur = trimmed.Ur
    Vr = trimmed.Vr

    (nt, nz) = size(w)
    A = Vector{Union{Missing, Float64}}(missing, nz)
    noise = Vector{Union{Missing, Float64}}(missing, nz)
    nrho  = zeros(Int, nz)
    R2_fit = Vector{Union{Missing, Float64}}(missing, nz)
    RMSE_fit = Vector{Union{Missing, Float64}}(missing, nz)
    nbins_fit = zeros(Int, nz)
    se_A_fit = Vector{Union{Missing, Float64}}(missing, nz)
    A_ci_lo = Vector{Union{Missing, Float64}}(missing, nz)
    A_ci_hi = Vector{Union{Missing, Float64}}(missing, nz)
    epsi_ci_lo = Vector{Union{Missing, Float64}}(missing, nz)
    epsi_ci_hi = Vector{Union{Missing, Float64}}(missing, nz)

    rhobin = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)
    D2bin  = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)
    rhovarbin = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)
    D2varbin = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)

    # Item 1 and 2: reuse pair indices and pass precomputed index vectors into displacements.
    level_cache = get_sf_level_cache(nt, nz)
    # Item 3: precompute scaled wind once per call.
    Udt = Ur * timestep
    Vdt = Vr * timestep

    for izo in 1:nz
        lc = level_cache[izo]
        rhoscale_max = max(Lscale_max, 5.0 * rng(izo))^(2/3)
        zm, dr2, dz2, D2, var_dr2, var_dz2 = displacements(
            lc.ci1, lc.ci2, lc.it1, lc.iz1, lc.it2, lc.iz2, 
            Udt, Vdt, pitch, roll, w, σ_pitch, σ_roll )
        rho = rhopair.(dr2, dz2) # r^2/3 * (1-(dz/dr)^2/4) for each pair
        # roughly exclude energy containing scales
        jj = map(r -> !ismissing(r) 
                 && isfinite(r) 
                 && (r <= rhoscale_max), rho)
        rho_var = propagate_rho_uncertainty.(
            dr2[jj], dz2[jj], var_dr2[jj], var_dz2[jj])

        # total number of valid realizations
        nrho[izo] = sum(jj .& .!ismissing.(D2) .& isfinite.(D2))

        (nbin_actual, rhobin_, D2inbin_, rhoinbin_, 
         D2varbin_, rhovarbin_) = equal_bin(
            rho[jj], D2[jj], rho_var; nbin_out_max=nbin_out_max)
        rhobin[1:nbin_actual, izo] .= rhoinbin_ # binned mean rho values not bin edge or centers rhobin_
        D2bin[1:nbin_actual, izo] .= D2inbin_
        rhovarbin[1:nbin_actual, izo] .= rhovarbin_
        D2varbin[1:nbin_actual, izo] .= D2varbin_
        
        # good values to condition matrices properly
        ii = (    .!ismissing.(rhobin[1:nbin_actual, izo])
               .& .!ismissing.(D2bin[1:nbin_actual, izo])
               .& .!ismissing.(rhovarbin[1:nbin_actual, izo])
               .&  ( rhovarbin[1:nbin_actual, izo] .> 0 )
               .& .!ismissing.(D2varbin[1:nbin_actual, izo])
               .&  ( D2varbin[1:nbin_actual, izo] .> 0 )
               )
        if sum(ii) >= 5
            if method == "tls"
                # total least squares fit of D2 = A * rho + noise, with errors in both variables
                try 
                    A_, noise_ = tls_slope_intercept(
                                        rhobin[1:nbin_actual, izo][ii], 
                                        D2bin[1:nbin_actual, izo][ii], 
                                        rhovarbin[1:nbin_actual, izo][ii], 
                                        D2varbin[1:nbin_actual, izo][ii],
                                        min_r_var=1e-2,
                                        min_d_var=1e-6 )
                    A[izo] = A_
                    noise[izo] = noise_
                catch e
                    nrmissing = sum(ismissing.(rhovarbin[1:nbin_actual, izo][ii])) # should be 0
                    ndmissing = sum(ismissing.(D2varbin[1:nbin_actual, izo][ii])) # should be 0
                    nrzero = sum((rhovarbin[1:nbin_actual, izo][ii].==0)) # should be 0
                    ndzero = sum((D2varbin[1:nbin_actual, izo][ii].==0)) # should be 0
                    @warn "TLS failed for level $izo with $nrmissing missing $nrzero zero rhovarbin and $ndmissing missing $ndzero zero D2varbin"
                end
            end
        end
        stats = fit_stats_onepass(rhoinbin_, D2inbin_)

        nbins_fit[izo] = stats.nbins
        # A[izo] = stats.A
        # noise[izo] = stats.noise
        R2_fit[izo] = stats.R2
        RMSE_fit[izo] = stats.RMSE
        se_A_fit[izo] = stats.se_A
        A_ci_lo[izo] = stats.A_ci_lo
        A_ci_hi[izo] = stats.A_ci_hi
        epsi_ci_lo[izo] = stats.epsi_ci_lo
        epsi_ci_hi[izo] = stats.epsi_ci_hi
    end
    return D2bin, rhobin, A, noise, R2_fit, RMSE_fit, nbins_fit, se_A_fit, A_ci_lo, A_ci_hi, epsi_ci_lo, epsi_ci_hi, nrho
end

# epsilon fitting using only small displacements
# defines functions small_displacements, small_bin, D2_rho_small_disp

"""
zm, dr2, dz2, D2, var_dr2, var_dz2, r_bin_edges = small_displacements(ci1,ci2,it1,iz1,it2,iz2,Udt,Vdt,pitch,roll,w; σ_pitch=σ_pitch, σ_roll=σ_roll, rangegate=rangegate)
Computes structure function D2 for
Displacements of sample pairs for one (vertical) level using precomputed pair-index vectors, 
including variance estimates.
"""
function small_displacements(ci1, ci2, it1, iz1, it2, iz2,
            Udt::Float64, Vdt::Float64, pitch, roll, w,
            σ_pitch::Real, σ_roll::Real; 
            rangegate=rangegate, Lmax=125.0 )

    # simplified to 1-level so iz1==iz2, dz=0
    rng(iz) = rangegate * (iz - 1 + 0.5) # range function

    # limit to small displacements in the inertial range
    # Lmax2 = Lmax^2
    dS = sqrt(Udt^2 + Vdt^2)   # scalar wind displacement for one time step
    itmax = floor(Int, Lmax/dS) # scalar maximum index
    r_bin_edges = dS * ((0:itmax) .+ 0.5) # bin by displacements later

    ii = @. abs(it2 - it1) <= itmax
    ni = count(ii)
    var_dr2 = Vector{Float64}(undef, ni)
    var_dz2 = Vector{Float64}(undef, ni)

    # view subset of ni valid indices
    let it1 = view(it1, ii), it2 = view(it2, ii), 
        iz1 = view(iz1, ii), iz2 = view(iz2, ii), 
        ci1 = view(ci1, ii), ci2 = view(ci2, ii)
        
        # displacement vectors of combinations of pairs
        X = @. Udt * (it2 - it1)
        Y = @. Vdt * (it2 - it1)
        zm = @. (rng(iz2)*cos(pitch[it2])*cos(roll[it2]) + rng(iz1)*cos(pitch[it1])*cos(roll[it1])) / 2
        dx = @. X + rng(iz2)*(-sin(pitch[it2])) - rng(iz1)*(-sin(pitch[it1]))
        dy = @. Y + rng(iz2)*cos(pitch[it2])*sin(roll[it2]) - rng(iz1)*cos(pitch[it1])*sin(roll[it1])
        dr2 = @. dx*dx + dy*dy
        dz2 = 0.0 # scalar!

        # structure function
        D2  = @. (w[ci2] - w[ci1])^2

        # Optimized loop using explicit 1:ni to bypass view-indexing overhead
        for i in 1:ni
            idx_t1 = it1[i]
            idx_t2 = it2[i]
            
            var_dr2[i], var_dz2[i] = DopplerTurbulence.displacement_variance_pitch_roll(
                rng(iz1[i]), rng(iz2[i]),
                pitch[idx_t1], pitch[idx_t2],
                roll[idx_t1], roll[idx_t2],
                σ_pitch, σ_roll               )
        end
        zm, dr2, dz2, D2, var_dr2, var_dz2, r_bin_edges # return these from function
    end
end

"""
nbin, ninbin, D2inbin, rhoinbin, D2varbin, rhovartot = small_bin(rho, D2, rho_var, dr2, r_bin_edges2) 
bin mean and variance of rho and D2 by dr2 within r_bin_edges2
Adds propagated instrument uncertainty to rhovarbin
"""
function small_bin(rho, D2, rho_var, dr2, r_bin_edges2)
    nbin = length(r_bin_edges2) - 1
    # rho_bin_edges = r_bin_edges2.^(2/3)
    D2inbin    = zeros(nbin)
    D2varbin   = zeros(nbin)
    rhoinbin   = zeros(nbin)
    rhovarbin  = zeros(nbin)
    rhovarmean = zeros(nbin)
    rhovartot  = zeros(nbin)
    ninbin     = zeros(Int, nbin)

    rho_var_finite = coalesce.(rho_var, 0.0)  # missing → 0 (no instrument error)
    ii = findall(.!ismissing.(rho) .& isfinite.(rho) .& .!ismissing.(D2) .& isfinite.(D2))
    if length(ii) >= 8
        for i in ii
            # index into bins
            idx = searchsortedlast(r_bin_edges2, dr2[i])
            if 1 <= idx <= nbin
                D2inbin[idx]    += D2[i]
                D2varbin[idx]   += D2[i]^2
                rhoinbin[idx]   += rho[i]
                rhovarbin[idx]  += rho[i]^2
                rhovarmean[idx] += rho_var_finite[i]
                ninbin[idx]     += 1
            end
        end
        D2inbin    ./= ninbin  # mean
        rhoinbin   ./= ninbin
        rhovarmean ./= ninbin
        D2varbin  = D2varbin  ./ ninbin - D2inbin.^2   # variance
        rhovarbin = rhovarbin ./ ninbin - rhoinbin.^2
        rhovartot = rhovarbin .+ rhovarmean

        return nbin, ninbin, D2inbin, rhoinbin, D2varbin, rhovartot
    else
        D2inbin    .= NaN
        rhoinbin   .= NaN
        D2varbin   .= NaN
        rhovartot  .= NaN
        ninbin     .= 0
        return nbin, ninbin, D2inbin, rhoinbin, D2varbin, rhovartot
    end
end

# ignore the smallest rho bin:
# dS = sqrt(Udt^2 + Vdt^2)
# Lmin = max(2*rangegate, 2*dS)
# r_bin_edges2 >= Lmin^2
# use total least squared to get
# regression slope of D2inbin vs. rhoinbin

"""
multiply D2 by inv_transfer_fcn(r2, dS2) to compensate for finite 
sample volume
"""
function inv_transfer_fcn(r2, dS2; rangegate2=rangegate^2)
    zi = (1.22/2.0)^2 # Gaussian window variance factor
    xi = 1/12 # rectangular window variance factor
    width2 = xi*dS2 + zi*rangegate2
    (1.0 + width2/r2)^(1/3)
end
"gaussian transfer function"
function transfer_fcn(r2, dS2, rangegate2=rangegate^2)
    1 - exp(-0.5*(r2/(dS2 + rangegate2))^(2/3))
end

"""
D2_rho_small_disp(w, pitch, roll, Ur, Vr, σ_pitch, σ_roll; nbin_out_max=25)
Replaces D2_rho_stare_qc, but also returns fit-quality metrics and confidence intervals.
No gating is applied; A and epsilon behavior remains unchanged.
"""
function D2_rho_small_disp(w, pitch, roll, Ur, Vr, 
                         σ_pitch::Real, σ_roll::Real;
                         nbin_out=100, method="tls", Lscale_max=3500.0)
    
    # paramters for good D2 for slope fits
    k0 = 3
    k1 = 8
    fac = 0.6

    trimmed = trim_structure_inputs(w, pitch, roll, Ur, Vr)
    nz_out = size(w, 2)
    out_A = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_noise = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_nrho = zeros(Int, nz_out)
    out_R2 = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_RMSE = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_nbins = zeros(Int, nz_out)
    out_se_A = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_A_lo = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_A_hi = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_epsi_lo = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_epsi_hi = Vector{Union{Missing, Float64}}(missing, nz_out)
    out_rho = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz_out)
    out_D2 = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz_out)

    if isnothing(trimmed)
        return out_D2, out_rho, out_A, out_noise, out_R2, out_RMSE, out_nbins, out_se_A, out_A_lo, out_A_hi, out_epsi_lo, out_epsi_hi, out_nrho
    end

    w = trimmed.w
    pitch = trimmed.pitch
    roll = trimmed.roll
    Ur = trimmed.Ur
    Vr = trimmed.Vr

    (nt, nz) = size(w)
    A = Vector{Union{Missing, Float64}}(missing, nz)
    noise = Vector{Union{Missing, Float64}}(missing, nz)
    nrho  = zeros(Int, nz)
    R2_fit = Vector{Union{Missing, Float64}}(missing, nz)
    RMSE_fit = Vector{Union{Missing, Float64}}(missing, nz)
    nbins_fit = zeros(Int, nz)
    se_A_fit = Vector{Union{Missing, Float64}}(missing, nz)
    A_ci_lo = Vector{Union{Missing, Float64}}(missing, nz)
    A_ci_hi = Vector{Union{Missing, Float64}}(missing, nz)
    epsi_ci_lo = Vector{Union{Missing, Float64}}(missing, nz)
    epsi_ci_hi = Vector{Union{Missing, Float64}}(missing, nz)

    rhobin = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)
    D2bin  = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)
    rhovarbin = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)
    D2varbin = Matrix{Union{Missing, Float64}}(missing, nbin_out, nz)

    # Item 1 and 2: reuse pair indices and pass precomputed index vectors into displacements.
    level_cache = get_sf_level_cache(nt, nz)
    # Item 3: precompute scaled wind once per call.
    Udt = Ur * timestep # vert vector
    Vdt = Vr * timestep
    dS2 = Udt.^2 .+ Vdt.^2 # vector

    for izo in 1:nz # loop vertical
        lc = level_cache[izo]
        rhoscale_max = max(Lscale_max, 5.0 * rng(izo))^(2/3) #superfluous

        zm, dr2, dz2, D2, var_dr2, var_dz2, r_bin_edges = small_displacements(
            lc.ci1, lc.ci2, lc.it1, lc.iz1, lc.it2, lc.iz2, 
            Udt[izo], Vdt[izo], pitch, roll, w, σ_pitch, σ_roll )

        rho = rhopair.(dr2, dz2) # r^2/3 * (1-(dz/dr)^2/4) for each pair
        # filter and coarsely exclude energy containing scales
        # (already displacements from small_displacements are stricter)
        jj = map(r -> !ismissing(r) 
                 && isfinite(r) 
                 && (r <= rhoscale_max), rho)
        zerodz = zeros(sum(jj)) # dz2 and var_dz2 are scalar zeros for small displacement version
        rho_var = propagate_rho_uncertainty.(
            dr2[jj], zerodz, var_dr2[jj], zerodz )
            # dr2[jj], dz2, var_dr2[jj], var_dz2 )

        # total number of valid realizations
        nrho[izo] = sum(jj .& .!ismissing.(D2) .& isfinite.(D2))

        (nbin_actual_, n_inbin, D2inbin_, rhoinbin_, D2varbin_, rhovarbin_) = small_bin(
            rho[jj], D2[jj], rho_var, dr2, r_bin_edges.^2 )
        # rhobin_ centers eliminated
        nbin_actual = min(nbin_actual_, nbin_out)
        # println("nbin_actual_=$(nbin_actual_), nbin_out=$(nbin_out), nbin_actual=$(nbin_actual)")
        rhobin[1:nbin_actual, izo] .= rhoinbin_[1:nbin_actual] # bin-mean rho values not bin edge or centers 
        # TXF
        D2bin[1:nbin_actual, izo] .= D2inbin_[1:nbin_actual] # .* inv_transfer_fcn.(rhoinbin_ .^3, 0.0) # dS2[izo])
        rhovarbin[1:nbin_actual, izo] .= rhovarbin_[1:nbin_actual]
        D2varbin[1:nbin_actual, izo] .= D2varbin_[1:nbin_actual]
        
        # bins for slope-intercept fitting
        ks = k0:min(k1, nbin_actual, round(Int, fac*argmax(D2bin[1:nbin_actual, izo])))

        # good values to condition matrices properly
        good(x) = !ismissing(x) && isfinite(x)
        ii = (    good.(rhobin[ks, izo])
               .& good.(D2bin[ks, izo])
               .& good.(rhovarbin[ks, izo]).& ( rhovarbin[ks, izo] .> 0 )
               .& good.(D2varbin[ks, izo] ).& ( D2varbin[ks, izo] .> 0 )
               )
        if sum(ii) >= 5
            if method == "tls"
                # total least squares fit of D2 = A * rho + noise, with errors in both variables
                try 
                    A_, noise_ = tls_slope_intercept(
                                        rhobin[ks, izo][ii], 
                                        D2bin[ks, izo][ii], 
                                        rhovarbin[ks, izo][ii], 
                                        D2varbin[ks, izo][ii],
                                        min_r_var=1e-2,
                                        min_d_var=1e-6 )
                    A[izo] = A_
                    noise[izo] = noise_
                catch e
                    nrmissing = sum(ismissing.(rhovarbin[ks, izo][ii])) # should be 0
                    ndmissing = sum(ismissing.(D2varbin[ks, izo][ii])) # should be 0
                    nrzero = sum((rhovarbin[ks, izo][ii].==0)) # should be 0
                    ndzero = sum((D2varbin[ks, izo][ii].==0)) # should be 0
                    @warn "TLS failed for level $izo with $nrmissing missing $nrzero zero rhovarbin and $ndmissing missing $ndzero zero D2varbin"
                end
            end
        end
        stats = fit_stats_onepass(rhoinbin_, D2inbin_)

        nbins_fit[izo] = stats.nbins
        # A[izo] = stats.A
        # noise[izo] = stats.noise
        R2_fit[izo] = stats.R2
        RMSE_fit[izo] = stats.RMSE
        se_A_fit[izo] = stats.se_A
        A_ci_lo[izo] = stats.A_ci_lo
        A_ci_hi[izo] = stats.A_ci_hi
        epsi_ci_lo[izo] = stats.epsi_ci_lo
        epsi_ci_hi[izo] = stats.epsi_ci_hi
    end
    return D2bin, rhobin, A, noise, R2_fit, RMSE_fit, nbins_fit, se_A_fit, A_ci_lo, A_ci_hi, epsi_ci_lo, epsi_ci_hi, nrho
end

# compute structure function slope
# parameters
k0 = 3
k1 = 8
# separations are spaced evenly in r, not rho
fac = 0.6 # r-factor from max D2 to avoid energy-containing range
# noisequantile = 0.2 # quantile of noise to subtract from D2

# functions
"""
noise from vertical mean intercept of D2 vs rho in inertial range k0:k1
not used
"""
function get_noise(rhobin, D2bin)
    noise = NaN .+ zeros(size(D2bin, 2))
    slope = NaN .+ zeros(size(D2bin, 2))
    for iz in iiz # loop over heights
        # good indices for linear fit at each height iiz
        kmx = argmax(D2bin[k0:k1+1, iz]) - 1
        ks = (k0:k1)[1:kmx]
        # linear fit D2 vs rho
        # p = polyfit(rhobin[ks, iz], D2bin[ks, iz], 1)
        slope[iz] = rhobin[ks, iz] \ D2bin[ks, iz]
        noise[iz] = mean(D2bin[ks, iz]) - slope[iz] * mean(rhobin[ks, iz])
    end
    noise, slope
end

"""
slope_ratio, slope_fit = get_slope(rho, D2; fac=fac)
slope D2/rho at each height from linear fit and simple ratio.
refines inertial range for finding slope.
remove noise from D2 first!
Since noise is not robustly determined, use slope_fit from linear regression.
"""
function get_slope(rho, D2; fac=fac)
    slopefit = NaN .+ zeros(size(D2, 2))
    slope = NaN .+ zeros(size(D2, 2))
    for iz in iiz # loop over heights
        kmx = round(Int, argmax(D2[:,iz]) * fac )
        ks = k0:min(k1, kmx) # range over which to calculate slope
        # fit slope from linear regression of D2 vs rho
        slopefit[iz] = rho[ks, iz] \ D2[ks, iz] # slope from linear regression
        # use simple ratio D2/rho
        # slope[iz] = median( D2[ks, iz] ./rho[ks, iz] )
    end
    slopefit
end

"get epsilon from slope A at each height"
get_epsilon(D2bin, rhobin) = epsilon.(get_slope( rhobin, D2bin ))

end # module DopplerTurbulence