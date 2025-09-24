"""
Julia routine that infers platform heave directly from the vertical stare by exploiting the fact that heave appears as a nearly range-uniform velocity shift. It is robust to rain contamination because it:

estimates a per-frame weighted robust location (weighted median) of the velocity profile using weights derived from backscatter (so aerosol-like gates weigh more),

uses that as an initial heave, subtracts it, builds a quick rain mask from the residuals + β anomalies,

recomputes weights excluding rain and re-estimates heave, and iterates to convergence.

No spectral width required. Uses only Base + Statistics.
"""

using Statistics

# -------------------------
# weighted median (ignores NaNs)
# -------------------------
"""
weighted_median(values::AbstractVector, weights::AbstractVector)

Return the weighted median of values with nonnegative weights.
NaNs in values or weights are ignored. If total weight is zero, returns NaN.
"""
function weighted_median(values::AbstractVector, weights::AbstractVector)
    @assert length(values) == length(weights)
    # collect valid pairs
    idx = findall(!((isnan∘identity).(values)) .& .!((isnan∘identity).(weights)) .& (weights .> 0))
    if isempty(idx)
        return NaN
    end
    vals = collect(values[idx])
    w    = collect(weights[idx])
    # sort by value
    perm = sortperm(vals)
    valss = vals[perm]
    ws = w[perm]
    total = sum(ws)
    half = total/2
    csum = 0.0
    for i in 1:length(ws)
        csum += ws[i]
        if csum >= half
            return valss[i]
        end
    end
    return valss[end]  # fallback
end

# -------------------------
# local vertical median & MAD for beta (used to compute β z-score)
# -------------------------
function local_vertical_stats(beta::AbstractMatrix; window::Int=11)
    R,T = size(beta)
    halfw = div(window,2)
    beta_med = fill(NaN, R, T)
    beta_mad = fill(NaN, R, T)
    for r in 1:R
        r0 = max(1, r-halfw); r1 = min(R, r+halfw)
        for t in 1:T
            seg = beta[r0:r1, t]
            segc = collect(skipmissing(seg))
            if isempty(segc)
                continue
            end
            m = median(segc)
            beta_med[r,t] = m
            beta_mad[r,t] = median(abs.(segc .- m))
            if beta_mad[r,t] == 0
                # small fallback to std
                beta_mad[r,t] = max(1e-12, std(segc))
            end
        end
    end
    return beta_med, beta_mad
end

# -------------------------
# Heave inference by iterative weighted-median & rain masking
# -------------------------
"""
infer_heave_from_stare(v::AbstractMatrix, beta::AbstractMatrix; kwargs...) -> (heave, rain_mask)

Inputs:
 - v[R,T] : Doppler vertical- stare velocities (m/s), NaN for missing gates.
 - beta[R,T]: backscatter (linear units), NaN for missing gates.

Returns:
 - heave[T] : inferred platform heave time series (m/s)
 - rain_mask[R,T] : boolean mask of rain-affected gates
Options (kwargs):
 - prefer_high_beta::Union{Nothing,Bool} : if nothing, auto-detect whether aerosol => high beta
 - init_beta_pctile = 50.0 : percentile used to form initial aerosol weight split
 - v_thresh = 0.6 : m/s threshold for fall-speed evidence
 - beta_z_thresh = 1.0 : z-score threshold for treating beta as anomalous
 - vert_min = 3 : minimum contiguous vertical length to accept rain
 - max_iter = 6 : max iterations
 - tol = 1e-3 : convergence tol on heave (m/s)
 - weight_power = 2.0 : how strongly to upweight aerosol-like gates (weights = (beta_rel)^power)
"""
function infer_heave_from_stare(v::AbstractMatrix, beta::AbstractMatrix; 
                                prefer_high_beta::Union{Nothing,Bool}=nothing,
                                init_beta_pctile=50.0,
                                v_thresh=0.6,
                                beta_z_thresh=1.0,
                                vert_min=3,
                                max_iter=6,
                                tol=1e-3,
                                weight_power=2.0,
                                beta_window=11)

    R,T = size(v)
    @assert size(beta) == (R,T)

    # compute vertical local beta median & mad for anomaly detection
    beta_med, beta_mad = local_vertical_stats(beta; window=beta_window)
    beta_z = (beta .- beta_med) ./ beta_mad

    # auto-detect whether aerosol corresponds to high or low beta
    if prefer_high_beta === nothing
        # choose gates with small temporal variance as aerosol candidates
        var_v = mapslices(x->var(collect(skipmissing(x))), v; dims=2)
        var_v = vec(var_v)
        nonnan_vars = var_v[.!isnan.(var_v)]
        if isempty(nonnan_vars)
            prefer_high_beta = true
        else
            thr = quantile(nonnan_vars, 0.3)   # low-variance gates
            cand = findall(x->(!isnan(x) && x <= thr), var_v)
            if isempty(cand)
                prefer_high_beta = true
            else
                # aggregate median beta for candidates vs global median
                cand_beta = collect(skipmissing(vec(beta[cand, :])))
                if isempty(cand_beta)
                    prefer_high_beta = true
                else
                    prefer_high_beta = median(cand_beta) >= median(collect(skipmissing(vec(beta))))
                end
            end
        end
    end

    # initial weights from beta at each frame: higher weight for aerosol-like gates
    # We'll form w_r,t = f( beta[r,t] relative to median at that range ), power-law
    # If prefer_high_beta true -> aerosol higher beta -> weight = max(0, beta - med)
    # Else -> aerosol lower beta -> weight = max(0, med - beta)
    base_weights = zeros(R,T)
    if prefer_high_beta
        base_weights = max.(0.0, beta .- beta_med)
    else
        base_weights = max.(0.0, beta_med .- beta)
    end
    # small floor so weight isn't exactly zero everywhere
    base_weights .= base_weights .+ 1e-6

    # normalize per-frame and apply power to sharpen
    weights = similar(base_weights)
    for t in 1:T
        w = base_weights[:,t]
        # clamp NaNs to zero
        w[isnan.(w)] .= 0.0
        s = sum(w)
        if s == 0
            # fallback: uniform weights for available gates
            w .= .!isnan.(v[:,t]) .* 1.0
            s = sum(w)
            if s == 0
                weights[:,t] .= 0.0
                continue
            end
        end
        w .= w ./ s
        w .= w .^ weight_power
        # renormalize
        sw = sum(w)
        if sw > 0
            w ./= sw
        end
        weights[:,t] = w
    end

    # iterative loop:
    heave = fill(NaN, T)
    prev_heave = fill(Inf, T)
    rain_mask = falses(R,T)

    for iter in 1:max_iter
        # compute weighted median per frame using current weights, ignoring gates flagged as rain
        for t in 1:T
            # mask gates currently identified as rain
            valid_idx = .!rain_mask[:,t] .& .!isnan.(v[:,t]) .& (weights[:,t] .> 0)
            if sum(valid_idx) == 0
                # fallback: use unweighted median of valid gates
                heave[t] = isempty(collect(skipmissing(v[:,t]))) ? NaN : median(collect(skipmissing(v[:,t])))
            else
                heave[t] = weighted_median(v[valid_idx, t], weights[valid_idx, t])
            end
        end

        # center velocities by heave and compute residuals
        vcorr = copy(v)
        for t in 1:T
            if !isnan(heave[t])
                vcorr[:,t] .-= heave[t]
            end
        end

        # quick rain scoring: evidence = (|vcorr| >= v_thresh) & (|beta_z| >= beta_z_thresh OR large local abs beta dev)
        candidate = falses(R,T)
        for r in 1:R, t in 1:T
            if isnan(vcorr[r,t]) || isnan(beta_z[r,t])
                continue
            end
            if abs(vcorr[r,t]) >= v_thresh && abs(beta_z[r,t]) >= beta_z_thresh
                candidate[r,t] = true
            end
            # also allow strong velocity alone if vertical coherence present (neighbor agreement)
        end

        # vertical coherence enforcement: require contiguous run length >= vert_min
        new_mask = falses(R,T)
        for t in 1:T
            r = 1
            while r <= R
                if candidate[r,t]
                    j = r
                    while j <= R && candidate[j,t]
                        j += 1
                    end
                    run = j - r
                    if run >= vert_min
                        new_mask[r:j-1, t] .= true
                    end
                    r = j
                else
                    r += 1
                end
            end
        end

        # narrow: also mark gates with |vcorr| >> v_thresh even if beta_z weak, but only if contiguous
        strongv = abs.(vcorr) .>= (1.5*v_thresh)
        for t in 1:T
            r = 1
            while r <= R
                if strongv[r,t]
                    j = r
                    while j <= R && strongv[j,t]
                        j += 1
                    end
                    run = j - r
                    if run >= vert_min
                        new_mask[r:j-1, t] .= true
                    end
                    r = j
                else
                    r += 1
                end
            end
        end

        # update rain_mask
        rain_mask .= new_mask

        # recompute weights suppressing rain gates strongly (set their weight to tiny epsilon)
        for t in 1:T
            w = base_weights[:,t]
            w[rain_mask[:,t]] .= 1e-8
            # clamp NaNs to zero
            w[isnan.(w)] .= 0.0
            s = sum(w)
            if s == 0
                w .= .!isnan.(v[:,t]) .* 1.0
                s = sum(w)
            end
            w .= w ./ s
            w .= w .^ weight_power
            sw = sum(w)
            if sw > 0
                w ./= sw
            end
            weights[:,t] = w
        end

        # check convergence of heave
        if iter > 1
            diffs = abs.(heave .- prev_heave)
            # ignore NaNs in max
            dmax = isempty(collect(skipmissing(diffs))) ? 0.0 : maximum(collect(skipmissing(diffs)))
            if dmax <= tol
                break
            end
        end
        prev_heave .= heave
    end

    return heave, rain_mask
end

