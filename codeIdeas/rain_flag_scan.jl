using Statistics
using LinearAlgebra
using ImageFiltering   # for median filter
using FFTW             # if you want later spectral ops
# If you don't have ImageFiltering, replace median with simple movmean or similar.

"""
rain_flag_scan(β, v; kwargs...) -> Bool array (mask) of same size as β and v.
Inputs:
 - β: backscatter array [range, azimuth]
 - v: Doppler velocity array [range, azimuth]
Both may contain NaN for missing gates.
Returns:
 - rainmask::BitArray{2} where true indicates rain-contaminated gate.

Options (kwargs):
 - β_low_pctile = 20.0    # percentile for "low backscatter" (since you said rain shows low aerosol β)
 - v_mag_thresh = 0.6     # m/s baseline fall-speed magnitude indicative of hydrometeors
 - az_consistency_frac = 0.25 # fraction of azimuth neighbors that must agree
 - min_vert_extent = 4    # minimum contiguous range gates to be considered a vertical band
 - median_kernel = (3,3)  # kernel for median smoothing
 - use_sigma_w = false    # optional spectral width array not included here
 - time_persistence = nothing # optional: supply previous mask or stack of masks
"""
function rain_flag_scan(β::AbstractMatrix, v::AbstractMatrix; 
                        β_low_pctile=20.0,
                        v_mag_thresh=0.6,
                        az_consistency_frac=0.25,
                        min_vert_extent=4,
                        median_kernel=(3,3),
                        temporal_history::Union{Nothing, Vector{BitArray{2}}}=nothing)
    R, A = size(β)
    # 1) Fill NaNs temporarily for percentile computations
    βflat = vec(β[.!isnan.(β)])
    if isempty(βflat)
        return falses(R,A)  # nothing to do
    end
    β_thresh = quantile(βflat, β_low_pctile/100)   # adaptive low-β threshold

    # 2) Pre-smooth β and v to suppress single-gate noise (median)
    # Use simple median filtering — ImageFiltering.medfilt works on arrays
    try
        βs = medfilt(β, median_kernel)
        vs  = medfilt(v, median_kernel)
    catch
        # fallback to original if medfilt not available
        βs = β
        vs  = v
    end

    # 3) Basic per-gate tests
    low_beta = map((b)->(!isnan(b) && b <= β_thresh), βs)
    strong_v  = map((vv)->(!isnan(vv) && abs(vv) >= v_mag_thresh), vs)

    # 4) Azimuthal consistency test:
    # For each range gate r, compute for azimuth a the fraction of neighbors (±N) that also show strong_v
    # We'll use a circular azimuth neighbor of width window_az
    window_az = max(3, Int(round(0.05*A)))  # ~5% of azimuths or at least 3
    halfw = div(window_az, 2)
    az_consistent = falses(R,A)
    for r in 1:R
        # precompute an array of Bool for this range
        row_strong = strong_v[r, :]
        for a in 1:A
            # circular window indices
            i1 = a - halfw
            inds = ((i1:i1+window_az-1) .- 1) .% A .+ 1
            frac = count(row_strong[inds]) / length(inds)
            az_consistent[r,a] = frac >= az_consistency_frac
        end
    end

    # 5) Vertical continuity test (vertical banding)
    # For each azimuth column a, find contiguous vertical runs where (low_beta & strong_v & az_consistent)
    candidate = low_beta .& strong_v .& az_consistent

    vert_band = falses(R,A)
    for a in 1:A
        r = 1
        while r <= R
            if candidate[r,a]
                # start run
                j = r
                while j <= R && candidate[j,a]
                    j += 1
                end
                runlen = j - r
                if runlen >= min_vert_extent
                    vert_band[r:j-1, a] .= true
                end
                r = j
            else
                r += 1
            end
        end
    end

    # 6) Small morphological cleaning: remove islands smaller than area_threshold
    # crude area removal: remove small connected patches (here we approximate with local sum filter)
    area_kernel = ones(Int, 5,5)   # 5x5 neighborhood
    # convert vert_band to Float to use convolution
    conv = imfilter(Float64.(vert_band), area_kernel; border="replicate")
    # area threshold: at least N cells (here 6)
    area_thresh = 6
    cleaned = vert_band .& (conv .>= area_thresh)

    # 7) Optional: temporal persistence — require mask to appear in multiple frames
    if temporal_history !== nothing
        # temporal_history is a vector of prior masks (BitArray) of same size
        n_required = max(1, ceil(Int, 0.5*(length(temporal_history)+1))) # require majority presence
        # form sum across history
        sum_hist = zeros(Int, R, A)
        for m in temporal_history
            sum_hist .+= Int.(m)
        end
        # include current cleaned as well
        sum_hist .+= Int.(cleaned)
        final_mask = sum_hist .>= n_required
    else
        final_mask = cleaned
    end

    return final_mask
end

