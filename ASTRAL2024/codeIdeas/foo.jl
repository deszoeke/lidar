using Statistics

# -------------------------
# Helper utilities
# -------------------------
# robust mean ignoring NaNs
robust_mean(x) = mean(skipmissing(x))
robust_median(x) = median(collect(skipmissing(x)))

# circular median along vector (no wrap needed for vertical stare)
function nanmedian(vec)
    v = collect(skipmissing(vec))
    isempty(v) && return NaN
    return median(v)
end

# simple 1D temporal smoothing (centered running median)
function running_median(x::AbstractVector, window::Int)
    n = length(x)
    out = fill(NaN, n)
    hw = div(window,2)
    for i in 1:n
        lo = max(1, i-hw); hi = min(n, i+hw)
        out[i] = nanmedian(x[lo:hi])
    end
    return out
end

# -------------------------
# 1) Estimate platform heave per time frame
#    Inputs:
#      v[range, time]  -> vertical stare Doppler velocity (NaN for bad gates)
#      beta[range, time] -> backscatter (NaN for bad gates)
#    Output:
#      heave[t] time series (m/s)
# Options:
#    - prefer_high_beta: if true assume aerosol produces high β; if false assume aerosol low β
#      We auto-detect by checking whether gates with low velocity variance have high or low β.
# -------------------------
function estimate_heave_ts(v::AbstractMatrix, beta::AbstractMatrix; 
                           prefer_high_beta::Union{Nothing,Bool}=nothing,
                           beta_pctile = 50.0,
                           variance_percentile = 30.0,
                           min_gates_for_heave=10,
                           temporal_smooth_win=3)

    R, T = size(v)
    @assert size(beta) == (R,T)

    # compute per-gate temporal variance of velocity (across time)
    var_v = mapslices(x->var(collect(skipmissing(x))), v; dims=2)  # R×1
    var_v = vec(var_v)

    # find low-variance gates (likely aerosol/clear-air where v ~ platform+air)
    varest_thresh = quantile(var_v[.!isnan.(var_v)], variance_percentile/100)
    candidate_gates = findall(x->(!isnan(x) && x <= varest_thresh), var_v)

    # decide whether aerosol corresponds to high or low beta
    if prefer_high_beta === nothing
        # compare median beta in candidate_gates vs full-median
        global_beta_med = quantile(vec(beta[.!isnan.(beta)]), 0.5)
        cand_beta = beta[candidate_gates, :]
        # compute median across time for candidate gates
        cand_beta_med = median(skipmissing(vec(cand_beta)))
        prefer_high_beta = cand_beta_med >= global_beta_med
    end

    # build per-frame heave estimate
    heave = fill(NaN, T)
    for t in 1:T
        # choose gates at time t that are temporal-low-var and also have appropriate beta
        idxs = candidate_gates
        if isempty(idxs)
            # fallback: use all gates with non-NaN
            vals = collect(skipmissing(v[:,t]))
            heave[t] = isempty(vals) ? NaN : median(vals)
            continue
        end
        betas_t = beta[idxs, t]
        vals_t  = v[idxs, t]
        # remove NaNs
        good = .!isnan.(betas_t) .& .!isnan.(vals_t)
        if sum(good) < min_gates_for_heave
            # widen selection: use any gate with low temporal variance and valid now
            vals_t2 = v[candidate_gates, t]
            vals = collect(skipmissing(vals_t2))
            heave[t] = isempty(vals) ? NaN : median(vals)
            continue
        end
        # select gates by beta relative to percentile at that time
        beta_t_all = beta[:,t]
        beta_pct = quantile(beta_t_all[.!isnan.(beta_t_all)], beta_pctile/100)
        if prefer_high_beta
            sel = (betas_t .>= beta_pct) .& good
        else
            sel = (betas_t .<= beta_pct) .& good
        end
        if sum(sel) < min_gates_for_heave
            # relax to use good
            sel = good
        end
        heave[t] = median(vals_t[sel])
    end

    # temporal smoothing to reduce spurious frame-to-frame jumps
    if temporal_smooth_win > 1
        heave = running_median(heave, temporal_smooth_win)
    end

    return heave, prefer_high_beta
end

# -------------------------
# 2) Rain scoring & masking for each (range, time) gate
#    score combines:
#      - beta anomaly (z-score relative to vertical local median)
#      - absolute heave-corrected velocity magnitude
#      - vertical coherence: number of adjacent range gates with similar sign/magnitude
#      - optional: spectral width (if provided)
#    Returns boolean rain_mask[range, time]
# -------------------------
function rain_mask_vertical(v::AbstractMatrix, beta::AbstractMatrix, heave::AbstractVector;
                            # thresholds and weights
                            v_mag_thresh=0.6,           # m/s - typical minimum fall speed
                            beta_z_thresh = 1.0,        # z-score threshold (positive or negative depending)
                            vert_coh_min = 3,           # min contiguous range gates
                            vert_coh_tol = 0.5,         # fraction agreement within window
                            spec_width = nothing,       # optional spectral width array [R,T]
                            spec_width_thresh = 0.5,    # optional
                            beta_window = 11,           # vertical window for local median
                            allow_beta_sign_flip = true # if optics ambiguous, allow low or high beta to indicate rain
                           )

    R, T = size(v)
    @assert size(beta) == (R,T)
    if spec_width !== nothing
        @assert size(spec_width) == (R,T)
    end

    # heave-corrected velocity
    vcorr = copy(v)
    for t in 1:T
        if !isnan(heave[t])
            vcorr[:,t] .-= heave[t]
        else
            # leave as is
        end
    end

    # compute local vertical median and MAD for beta
    halfw = div(beta_window,2)
    beta_med = Array{Float64}(undef, R, T)
    beta_mad = Array{Float64}(undef, R, T)
    for r in 1:R
        r0 = max(1, r-halfw); r1 = min(R, r+halfw)
        for t in 1:T
            seg = beta[r0:r1, t]
            segc = collect(skipmissing(seg))
            if isempty(segc)
                beta_med[r,t] = NaN; beta_mad[r,t] = NaN
            else
                m = median(segc)
                beta_med[r,t] = m
                beta_mad[r,t] = median(abs.(segc .- m))
                if beta_mad[r,t] == 0
                    beta_mad[r,t] = maximum([1e-12, std(segc)]) # avoid zero
                end
            end
        end
    end

    # z-score for beta: (beta - local_median)/mad
    beta_z = (beta .- beta_med) ./ beta_mad

    # vertical coherence: for each gate, fraction of neighbors within ±N that have abs(vcorr) >= v_mag_thresh and sign same
    vert_window = 2  # check ±2 gates (adjust with range gate spacing)
    vert_coh_frac = zeros(Float64, R, T)
    for r in 1:R
        r0 = max(1, r-vert_window); r1 = min(R, r+vert_window)
        for t in 1:T
            neigh = vcorr[r0:r1, t]
            me = vcorr[r,t]
            if isnan(me)
                vert_coh_frac[r,t] = 0.0
                continue
            end
            valid = .!isnan.(neigh)
            if sum(valid) == 0
                vert_coh_frac[r,t] = 0.0
                continue
            end
            # count neighbors with magnitude >= thresh and same sign as center
            same = abs.(neigh[valid]) .>= v_mag_thresh
            if !isnan(me)
                same .&= sign.(neigh[valid]) .== sign(me)
            end
            vert_coh_frac[r,t] = sum(same) / sum(valid)
        end
    end

    # assemble score: heavy weighting to v magnitude and vertical coherence; beta z sign handled below
    score = zeros(Float64, R, T)
    # normalized contributions
    vscore = clamp.(abs.(vcorr) ./ max(v_mag_thresh, 1e-6), 0.0, 4.0)  # values ~1 when at threshold
    cohscore = vert_coh_frac ./ 1.0
    score .= 0.6 .* vscore .+ 0.3 .* cohscore

    # beta evidence: allow either sign depending on local behavior. If allow_beta_sign_flip, treat large |z| as evidence.
    if allow_beta_sign_flip
        beta_evidence = clamp.(abs.(beta_z) ./ max(beta_z_thresh, 1e-6), 0.0, 2.0)
        score .+= 0.1 .* beta_evidence
    else
        # assume a particular sign indicates rain (user can change)
        # if positive z indicates rain
        beta_evidence = clamp.(beta_z ./ max(beta_z_thresh, 1e-6), 0.0, 2.0)
        score .+= 0.1 .* beta_evidence
    end

    # spec width optional
    if spec_width !== nothing
        spec_evi = clamp.(spec_width ./ max(spec_width_thresh, 1e-6), 0.0, 2.0)
        score .+= 0.05 .* spec_evi
    end

    # threshold score to produce mask; require vertical extent check
    raw_mask = score .>= 0.8   # tune this; 0.8 is a conservative starting point

    # require minimum vertical contiguous length vert_coh_min
    final_mask = falses(R, T)
    for t in 1:T
        r = 1
        while r <= R
            if raw_mask[r,t]
                j = r
                while j <= R && raw_mask[j,t]
                    j += 1
                end
                runlen = j - r
                if runlen >= vert_coh_min
                    final_mask[r:j-1, t] .= true
                end
                r = j
            else
                r += 1
            end
        end
    end

    return final_mask, score, vcorr
end

# -------------------------
# 3) Estimate rain fall speed per contiguous layer and subtract (optional)
#    For each time t and contiguous vertical block flagged as rain, compute robust mean/mode of vcorr
#    and subtract it from that block to recover ambient air vertical velocity.
# -------------------------
function subtract_rain_layers!(vcorr::AbstractMatrix, rain_mask::BitArray{2}; method=:median, min_count=3)
    R, T = size(vcorr)
    v_corrected = copy(vcorr)
    rain_layers = []  # collect info tuples (t, rstart, rend, vrain)
    for t in 1:T
        r = 1
        while r <= R
            if rain_mask[r,t]
                j = r
                while j <= R && rain_mask[j,t]
                    j += 1
                end
                runlen = j - r
                vals = collect(skipmissing(vcorr[r:j-1, t]))
                if length(vals) >= min_count
                    vrain = method == :median ? median(vals) : mean(vals)
                    # subtract vrain from that block
                    v_corrected[r:j-1, t] .-= vrain
                    push!(rain_layers, (t, r, j-1, vrain))
                end
                r = j
            else
                r += 1
            end
        end
    end
    return v_corrected, rain_layers
end

# -------------------------
# 4) Iterative wrapper: estimate heave, mask rain, recompute heave from clean gates, until converged
# -------------------------
function iterative_despike_rain(v::AbstractMatrix, beta::AbstractMatrix; niter=4, verbose=true)
    R, T = size(v)
    # initial heave estimate
    heave, prefer_high_beta = estimate_heave_ts(v, beta)
    for iter in 1:niter
        if verbose; println("Iteration $iter: estimating rain mask..."); end
        mask, score, vcorr = rain_mask_vertical(v, beta, heave)
        # mask rain and compute new heave from non-rain gates
        v_masked = copy(v)
        for t in 1:T
            v_masked[mask[:,t], t] .= NaN
        end
        new_heave, _ = estimate_heave_ts(v_masked, beta)
        # check convergence
        if all(isnan.(heave) .& isnan.(new_heave))
            heave = new_heave
            break
        end
        change = maximum(abs.(skipmissing(heave .- new_heave)))
        if verbose; println(" max heave change = $(isnan(change) ? 0.0 : change) m/s"); end
        heave = new_heave
    end

    # final mask with final heave
    final_mask, final_score, final_vcorr = rain_mask_vertical(v, beta, heave)
    # subtract rain falls
    v_corrected, layers = subtract_rain_layers!(final_vcorr, final_mask)
    return (heave=heave, rain_mask=final_mask, rain_score=final_score, v_corr=v_corrected, rain_layers=layers, prefer_high_beta=prefer_high_beta)
end

