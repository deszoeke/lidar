module LidarVNSync

using Dates
using Statistics
using StatsBase
using Rotations
using Interpolations
using DSP
using FFTW
using NCDatasets
using JLD2
using Printf
using PyPlot

include("./read_lidar.jl")
using .read_lidar
using .read_lidar.stare
using .read_vecnav: read_vecnav_dict
using .chunks: read_stare_chunk, read_stare_time, fit_offset

include("./timing_lidar.jl")
using .timing_lidar

export NOISE_THR, BANDPASS_PERIOD, BANDPASS_HZ, TIMESTEP, RANGEGATE
export init_periodic_beams, load_lidar_indices_and_files, read_streamlinexr_stare!
export init_stream_state, extract_sync_window, coarse_and_fine_lag, run_sequential_offsets, refine_offset_20hz
export setup_sync_context, process_sync_data, save_sync_netcdf
export diagnostic_single_window, diagnostic_offsets
export read_synced_motion, motion_correct_stare_velocity, read_and_motion_correct_stare

const NOISE_THR = 1.03
const BANDPASS_PERIOD = (5.0, 20.0)
const BANDPASS_HZ = (1 / BANDPASS_PERIOD[2], 1 / BANDPASS_PERIOD[1])
const TIMESTEP = 1.02
const RANGEGATE = 24.0

pd = permutedims

isgoodnum(x) = !ismissing(x) && isfinite(x)

function diag_array(name, A)
    v = vec(A)
    @printf("%-12s size=%-14s missing=%6d finite=%6d nan=%6d\n",
        name,
        string(size(A)),
        count(ismissing, v),
        count(isgoodnum, v),
        count(x -> !ismissing(x) && isnan(x), v),
    )
end

function valid_chunk_mdv(intensity, dopplervel; thr=NOISE_THR)
    valid = (.!ismissing.(dopplervel)) .& isfinite.(coalesce.(dopplervel, NaN)) .& (.!ismissing.(intensity)) .& (coalesce.(intensity, -Inf) .> thr)
    goodlevels = findall(vec(all(valid, dims=1)))
    nt = size(dopplervel, 1)
    if isempty(goodlevels)
        return fill(NaN, nt), goodlevels
    end
    chunkmedian = median(skipmissing(vec(dopplervel[:, goodlevels])))
    timeslicemedian = vec(mapslices(x -> median(skipmissing(x)), dopplervel[:, goodlevels], dims=2))
    return timeslicemedian .- chunkmedian, goodlevels
end

_wrap(x::Integer, n::Int) = mod1(x, n)
_wrap(x::AbstractRange{<:Integer}, n::Int) = mod1.(collect(x), n)
_wrap(x::AbstractArray{<:Integer}, n::Int) = mod1.(x, n)
_wrap(::Colon, ::Int) = (:)

struct PeriodicVector{T}
    data::Vector{T}
end

Base.size(p::PeriodicVector) = size(p.data)
Base.length(p::PeriodicVector) = length(p.data)
Base.getindex(p::PeriodicVector, i) = getindex(p.data, _wrap(i, length(p.data)))
Base.setindex!(p::PeriodicVector, v, i) = setindex!(p.data, v, _wrap(i, length(p.data)))
Base.setindex!(p::PeriodicVector, v::AbstractVector, i::Union{AbstractVector, AbstractRange}) = (p.data[_wrap(i, length(p.data))] .= v)

struct PeriodicMatrix{T}
    data::Matrix{T}
end

Base.size(p::PeriodicMatrix) = size(p.data)
Base.getindex(p::PeriodicMatrix, i, j...) = getindex(p.data, _wrap(i, size(p.data, 1)), j...)
Base.setindex!(p::PeriodicMatrix, v, i, j...) = setindex!(p.data, v, _wrap(i, size(p.data, 1)), j...)
Base.setindex!(p::PeriodicMatrix, v::AbstractArray, i::Union{AbstractVector, AbstractRange}, j...) = (p.data[_wrap(i, size(p.data, 1)), j...] .= v)

function init_periodic_beams(nbeams, ngates)
    PeriodicVectorMissing(n) = PeriodicVector(Vector{Union{Float32, Missing}}(fill(missing, n)))
    PeriodicMatrixMissing(nrows, ncols) = PeriodicMatrix(Matrix{Union{Float32, Missing}}(fill(missing, nrows, ncols)))
    VectorMissing(n) = Vector{Union{Float32, Missing}}(fill(missing, n))
    Dict(
        :time => PeriodicVectorMissing(nbeams),
        :azimuth => PeriodicVectorMissing(nbeams),
        :elevangle => PeriodicVectorMissing(nbeams),
        :pitch => PeriodicVectorMissing(nbeams),
        :roll => PeriodicVectorMissing(nbeams),
        :height => VectorMissing(ngates),
        :dopplervel => PeriodicMatrixMissing(nbeams, ngates),
        :intensity => PeriodicMatrixMissing(nbeams, ngates),
        :beta => PeriodicMatrixMissing(nbeams, ngates),
    )
end

function load_lidar_indices_and_files(; lidarstemdir="./data")
    if isfile("lidar_dt.jld2")
        lidardt = load("lidar_dt.jld2")
    else
        error("Expected lidar_dt.jld2 to exist. Generate it with save_lidar_dt.jl before using this workbench.")
    end
    if isfile("file_beam_inds.jld2")
        fileinds = load("file_beam_inds.jld2")
    else
        error("Expected file_beam_inds.jld2 to exist. Generate it from the production workflow before using this workbench.")
    end
    starefiles = filter(startswith("Stare_116_"), readdir(joinpath(lidarstemdir, "all")))
    ff = joinpath.(lidarstemdir, "all", starefiles)
    return (;
        LidarDt=lidardt,
        FileInds=fileinds,
        ff,
        ists=lidardt["ist"],
        iens=lidardt["ien"],
        dtime=lidardt["dtime"],
        bigind_file_starts=fileinds["bigind_file_start"],
        bigind_file_ends=fileinds["bigind_file_end"],
    )
end

function read_streamlinexr_stare!(file_path, h, beams, bb, nheaderlines=17; startat=1, endat=0)
    nz = size(beams[:height][:], 1)
    nlines = h[:nlines]
    ngates = h[:ngates]
    nbeamsmax = round(Int, (nlines - nheaderlines) / (1 + ngates))
    endat = endat == 0 ? nbeamsmax : mod(endat - 1, nbeamsmax) + 1
    nbeams = min(endat - startat + 1, nbeamsmax)

    beam_timeangles = zeros(Float64, nbeams, 5)
    beam_velrad = zeros(Float64, nbeams, ngates, 4)

    open(file_path) do file
        for _ in 1:nheaderlines
            readline(file)
        end
        for _ in 1:((1 + ngates) * (startat - 1))
            readline(file)
        end
        for ibeam in 1:nbeams
            beam_timeangles[ibeam, :] .= parse.(Float64, split(readline(file)))
            for igate in 1:ngates
                beam_velrad[ibeam, igate, :] .= parse.(Float64, split(readline(file)))
            end
        end
    end

    setindex!(beams[:time], beam_timeangles[:, 1], bb)
    setindex!(beams[:azimuth], beam_timeangles[:, 2], bb)
    setindex!(beams[:elevangle], beam_timeangles[:, 3], bb)
    setindex!(beams[:pitch], beam_timeangles[:, 4], bb)
    setindex!(beams[:roll], beam_timeangles[:, 5], bb)
    beams[:height][1:nz] .= (beam_velrad[1, 1:nz, 1] .+ 0.5) .* h[:gatelength]
    setindex!(beams[:dopplervel], beam_velrad[:, 1:nz, 2], bb, 1:nz)
    setindex!(beams[:intensity], beam_velrad[:, 1:nz, 3], bb, 1:nz)
    setindex!(beams[:beta], beam_velrad[:, 1:nz, 4], bb, 1:nz)
    nothing
end

function setup_sync_context(; lidarstemdir="./data", uv_path=joinpath("data/netcdf", "ekamsat_lidar_uv_20240428-20240613.nc"), nx=4000, nz=80)
    Env = load_lidar_indices_and_files(; lidarstemdir=lidarstemdir)
    Vn = read_vecnav_dict()
    UV = NCDataset(uv_path)

    dtime_st = Env.dtime[Env.ists]
    dtime_en = Env.dtime[Env.iens]
    icvn = findfirst(dtime_st .>= Vn[:vndt][1]):findlast(dtime_en .<= Vn[:vndt][end])
    beams = init_periodic_beams(nx, nz)

    return (; Env, Vn, UV, beams, icvn, nz, nx)
end

function init_stream_state()
    Dict(
        :ifile_loaded => 0,
        :bigind_file_end => 0,
    )
end

function ensure_chunk_loaded!(beams, env, state, ic)
    ist = env.ists[ic]
    ien = env.iens[ic]
    ifile_needed = findlast(env.bigind_file_starts .<= ist)
    isnothing(ifile_needed) && error("No file contains ist=$(ist)")

    if state[:ifile_loaded] == 0
        state[:ifile_loaded] = ifile_needed - 1
        state[:bigind_file_end] = 0
    end

    while (state[:ifile_loaded] < ifile_needed) || (ien > state[:bigind_file_end])
        state[:ifile_loaded] += 1
        state[:ifile_loaded] > length(env.ff) && error("Ran out of lidar files while loading chunk ic=$(ic)")

        ifile = state[:ifile_loaded]
        bb = env.bigind_file_starts[ifile]:env.bigind_file_ends[ifile]
        h = read_lidar.read_streamlinexr_head(env.ff[ifile])
        read_streamlinexr_stare!(env.ff[ifile], h, beams, bb)
        state[:bigind_file_end] = env.bigind_file_ends[ifile]
    end

    return (; ifile=ifile_needed, ist, ien)
end

function chunk_lidar_datetimes(dt_chunk, beams, ist, ien)
    stare_dt_raw = @. DateTime(Date(dt_chunk)) + Millisecond(round(Int64, beams[:time][ist:ien] * 3_600_000))
    lidar_clock_fast_by = Millisecond(round(Int64, 1_000 * fit_offset(stare_dt_raw[1])))
    stare_dt = stare_dt_raw .- lidar_clock_fast_by
    return stare_dt_raw, stare_dt, lidar_clock_fast_by
end

# UV is threaded through because read_stare_chunk is the shared dissipation-era reader
# and this extracted window will later feed motion correction. The sync offset fit below
# depends only on mdv and vn2, not on Ur/Vr.
function extract_sync_window(beams, env, state, Vn, UV, ic; ntop=80)
    info = ensure_chunk_loaded!(beams, env, state, ic)
    dt_arg = info.ist <= length(env.dtime) ? env.dtime[info.ist] : env.dtime[end]
    dopplervel, pitch, roll, vn0, vn1, vn2, Ur, Vr, mdv_builtin = read_stare_chunk(dt_arg, beams, Vn, UV, info.ist, info.ien, ntop)
    stare_dt_raw, stare_dt, lidar_clock_fast_by = chunk_lidar_datetimes(env.dtime[info.ist], beams, info.ist, info.ien)
    intensity = beams[:intensity][info.ist:info.ien, 1:ntop]
    beta = beams[:beta][info.ist:info.ien, 1:ntop]
    mdv, goodlevels = valid_chunk_mdv(intensity, dopplervel)
    return (; info..., stare_dt_raw, stare_dt, lidar_clock_fast_by, dopplervel, intensity, beta, pitch, roll, vn0, vn1, vn2, Ur, Vr, mdv, mdv_builtin, goodlevels)
end

function finite_overlap(x, y)
    n = min(length(x), length(y))
    xv = x[1:n]
    yv = y[1:n]
    keep = isfinite.(xv) .& isfinite.(yv)
    return xv[keep], yv[keep], keep
end

function detrend_center(x)
    y = Float64.(x)
    y .-= mean(y)
    y
end

function shift_signal_linear(x, lag_seconds; dt=TIMESTEP, fill_value=NaN)
    y = Float64.(x)
    n = length(y)
    q = collect(1:n) .- lag_seconds / dt
    out = fill(fill_value, n)
    for i in eachindex(q)
        qi = q[i]
        if 1 <= qi <= n
            lo = floor(Int, qi)
            hi = ceil(Int, qi)
            if lo == hi
                out[i] = y[lo]
            else
                w = qi - lo
                out[i] = (1 - w) * y[lo] + w * y[hi]
            end
        end
    end
    out
end

real_fft_frequencies(n, fs) = collect(0:div(n, 2)) .* (fs / n)

function cosine_edge_mask(freq, flo, fhi; frac=0.15)
    mask = zeros(Float64, length(freq))
    width = max((fhi - flo) * frac, eps())
    lo1 = max(0.0, flo - width)
    hi2 = fhi + width
    for i in eachindex(freq)
        f = freq[i]
        if flo <= f <= fhi
            mask[i] = 1.0
        elseif lo1 <= f < flo
            mask[i] = 0.5 * (1 - cos(pi * (f - lo1) / max(flo - lo1, eps())))
        elseif fhi < f <= hi2
            mask[i] = 0.5 * (1 + cos(pi * (f - fhi) / max(hi2 - fhi, eps())))
        end
    end
    mask
end

function fft_bandpass(x; dt=TIMESTEP, flo=BANDPASS_HZ[1], fhi=BANDPASS_HZ[2], taper_frac=0.1)
    n = length(x)
    n < 8 && return fill(NaN, n)
    xt = detrend_center(x)
    taper = DSP.Windows.tukey(n, taper_frac)
    xf = rfft(xt .* taper)
    freq = real_fft_frequencies(n, 1 / dt)
    mask = cosine_edge_mask(freq, flo, fhi)
    irfft(xf .* mask, n)
end

function analytic_envelope_fft(x)
    n = length(x)
    X = fft(Float64.(x))
    h = zeros(Float64, n)
    if iseven(n)
        h[1] = 1
        h[div(n, 2) + 1] = 1
        h[2:div(n, 2)] .= 2
    else
        h[1] = 1
        h[2:div(n + 1, 2)] .= 2
    end
    abs.(ifft(X .* h))
end

function fft_xcorr_lag(x, y; dt=TIMESTEP, maxlag_seconds=120.0, center_seconds=0.0)
    xv, yv, _ = finite_overlap(Float64.(x), Float64.(y))
    n = min(length(xv), length(yv))
    n < 8 && return (lag_seconds=NaN, peak=NaN, peak_norm=NaN, lags=Float64[], corr=Float64[])

    xuse = detrend_center(xv[1:n])
    yuse = detrend_center(yv[1:n])
    nfft = nextpow(2, 2 * n - 1)
    xpad = vcat(xuse, zeros(nfft - n))
    ypad = vcat(yuse, zeros(nfft - n))
    cc = ifft(fft(xpad) .* conj(fft(ypad)))
    corr = real(vcat(cc[end - (n - 2):end], cc[1:n]))
    lags = collect(-(n - 1):(n - 1)) .* dt

    keep = abs.(lags .- center_seconds) .<= maxlag_seconds
    corrw = corr[keep]
    lagsw = lags[keep]
    isempty(corrw) && return (lag_seconds=NaN, peak=NaN, peak_norm=NaN, lags=Float64[], corr=Float64[])

    peak, idx = findmax(corrw)
    denom = sqrt(sum(abs2, xuse) * sum(abs2, yuse))
    peak_norm = denom > 0 ? peak / denom : NaN
    return (lag_seconds=lagsw[idx], peak=peak, peak_norm=peak_norm, lags=lagsw, corr=corrw)
end

function iterative_coarse_lag(mdv_bp, vn2_bp; dt=TIMESTEP, prior_seconds=0.0, max_passes=4, search_seconds=12.0, tol_seconds=0.25, min_improve=1e-3)
    mdv_env = analytic_envelope_fft(mdv_bp)
    offset = prior_seconds
    history = NamedTuple[]
    prev_score = -Inf

    for pass in 1:max_passes
        vn2_bp_shifted = shift_signal_linear(vn2_bp, offset; dt=dt)
        vn2_env_shifted = analytic_envelope_fft(vn2_bp_shifted)
        coarse = fft_xcorr_lag(mdv_env, vn2_env_shifted; dt=dt, center_seconds=0.0, maxlag_seconds=search_seconds)
        improvement = pass == 1 || !isfinite(prev_score) ? Inf : coarse.peak_norm - prev_score
        push!(history, (; pass, offset_in=offset, residual=coarse.lag_seconds, peak=coarse.peak, peak_norm=coarse.peak_norm, improvement))

        if !isfinite(coarse.lag_seconds)
            break
        end

        offset += coarse.lag_seconds
        prev_score = coarse.peak_norm

        if abs(coarse.lag_seconds) <= tol_seconds
            break
        end
        if pass > 1 && isfinite(improvement) && improvement < min_improve
            break
        end
    end

    vn2_bp_coarse = shift_signal_linear(vn2_bp, offset; dt=dt)
    vn2_env_coarse = analytic_envelope_fft(vn2_bp_coarse)
    return (; coarse_offset=offset, mdv_env, vn2_bp_coarse, vn2_env_coarse, history)
end

robust_center(x) = isempty(x) ? NaN : median(x)

function prior_from_history(history, start_dt; recent_window=Minute(10))
    if isempty(history)
        return (; prior_seconds=0.0, previous_offset=NaN, recent_offset=NaN, parent_offset=NaN)
    end

    previous_offset = history[end].final_offset
    recent_offsets = [h.final_offset for h in history if start_dt - h.start_dt <= recent_window]
    parent_start = floor(start_dt, Hour)
    parent_offsets = [h.final_offset for h in history if floor(h.start_dt, Hour) == parent_start]

    recent_offset = robust_center(recent_offsets)
    parent_offset = robust_center(parent_offsets)

    values = Float64[]
    weights = Float64[]
    if isfinite(previous_offset)
        push!(values, previous_offset); push!(weights, 0.55)
    end
    if isfinite(recent_offset)
        push!(values, recent_offset); push!(weights, 0.30)
    end
    if isfinite(parent_offset)
        push!(values, parent_offset); push!(weights, 0.15)
    end

    prior_seconds = isempty(values) ? 0.0 : sum(values .* weights) / sum(weights)
    return (; prior_seconds, previous_offset, recent_offset, parent_offset)
end

function coarse_and_fine_lag(mdv, vn2; dt=TIMESTEP, prior_seconds=0.0, coarse_search_seconds=12.0, max_passes=4, fine_search_seconds=4.0)
    mdv_bp = fft_bandpass(mdv; dt=dt)
    vn2_bp = fft_bandpass(vn2; dt=dt)

    coarse = iterative_coarse_lag(mdv_bp, vn2_bp; dt=dt, prior_seconds=prior_seconds, max_passes=max_passes, search_seconds=coarse_search_seconds)
    fine = fft_xcorr_lag(mdv_bp, coarse.vn2_bp_coarse; dt=dt, center_seconds=0.0, maxlag_seconds=fine_search_seconds)
    final_offset = coarse.coarse_offset + fine.lag_seconds
    vn2_bp_final = shift_signal_linear(vn2_bp, final_offset; dt=dt)
    vn2_env_final = analytic_envelope_fft(vn2_bp_final)

    return (; prior_seconds, mdv_bp, vn2_bp, mdv_env=coarse.mdv_env, vn2_env=coarse.vn2_env_coarse, coarse_history=coarse.history, coarse_offset=coarse.coarse_offset, fine, final_offset, vn2_bp_coarse=coarse.vn2_bp_coarse, vn2_bp_final, vn2_env_final)
end

function is_jump_candidate(prior_seconds, final_seconds, fine_peak_norm; jump_threshold_seconds=0.8, min_fine_peak=0.08)
    isfinite(prior_seconds) && isfinite(final_seconds) && isfinite(fine_peak_norm) &&
    abs(final_seconds - prior_seconds) >= jump_threshold_seconds && fine_peak_norm >= min_fine_peak
end

function backward_jump_robustness(beams, env, state, Vn, UV, history, post_offset; ntop=80, n_back=2)
    isempty(history) && return NaN
    ntest = min(n_back, length(history))
    deltas = Float64[]
    for k in 1:ntest
        h = history[end - k + 1]
        win_prev = extract_sync_window(beams, env, state, Vn, UV, h.ic; ntop=ntop)
        back_sync = coarse_and_fine_lag(win_prev.mdv, win_prev.vn2; prior_seconds=post_offset)
        if isfinite(back_sync.final_offset) && isfinite(h.final_offset)
            push!(deltas, abs(back_sync.final_offset - h.final_offset))
        end
    end
    isempty(deltas) ? NaN : median(deltas)
end

function run_sequential_offsets(beams, env, Vn, UV, ic_list; ntop=80, jump_threshold_seconds=0.8, min_fine_peak=0.08, backward_windows=2, backward_tol=0.35)
    state = init_stream_state()
    history = NamedTuple[]
    for ic in ic_list
        start_dt = env.dtime[env.ists[ic]]
        prior = prior_from_history(history, start_dt)
        win = extract_sync_window(beams, env, state, Vn, UV, ic; ntop=ntop)
        sync = coarse_and_fine_lag(win.mdv, win.vn2; prior_seconds=prior.prior_seconds)

        jump_candidate = is_jump_candidate(prior.prior_seconds, sync.final_offset, sync.fine.peak_norm; jump_threshold_seconds=jump_threshold_seconds, min_fine_peak=min_fine_peak)

        jump_backward_metric = NaN
        jump_robust = false
        jump_accepted = false

        if jump_candidate && !isempty(history)
            post_sync = coarse_and_fine_lag(win.mdv, win.vn2; prior_seconds=sync.final_offset)
            post_offset = post_sync.final_offset
            jump_backward_metric = backward_jump_robustness(beams, env, state, Vn, UV, history, post_offset; ntop=ntop, n_back=backward_windows)
            jump_robust = isfinite(jump_backward_metric) && jump_backward_metric <= backward_tol
            if jump_robust
                sync = post_sync
                jump_accepted = true
            else
                sync = coarse_and_fine_lag(win.mdv, win.vn2; prior_seconds=prior.prior_seconds, fine_search_seconds=2.0)
            end
        end

        push!(history, (; ic, start_dt=win.stare_dt[1], end_dt=win.stare_dt[end], parent_start=floor(win.stare_dt[1], Hour), prior_offset=prior.prior_seconds, previous_offset=prior.previous_offset, recent_offset=prior.recent_offset, parent_offset=prior.parent_offset, coarse_offset=sync.coarse_offset, fine_residual=sync.fine.lag_seconds, final_offset=sync.final_offset, coarse_peak=isempty(sync.coarse_history) ? NaN : sync.coarse_history[end].peak_norm, fine_peak=sync.fine.peak_norm, goodlevels=length(win.goodlevels), jump_candidate, jump_robust, jump_accepted, jump_backward_metric))
    end
    history
end

function vn_subset_20hz(stare_dt, Vn; pad=Second(2))
    ind = findall(stare_dt[1] - pad .<= Vn[:vndt] .<= stare_dt[end] + pad)
    vndt = Vn[:vndt][ind]
    vn2 = Float64.(Vn[:VelNED2][ind])
    return ind, vndt, vn2
end

function dt_seconds(dtv::Vector{DateTime})
    length(dtv) < 2 && return NaN
    d = Float64.(Dates.value.(Millisecond.(diff(dtv)))) ./ 1000
    d = d[isfinite.(d) .& (d .> 0)]
    isempty(d) ? NaN : median(d)
end

function nan_safe_moving_average(x, n)
    out = fill(NaN, length(x))
    half = n ÷ 2
    for i in eachindex(x)
        lo = max(1, i - half)
        hi = min(length(x), i + half)
        s = 0.0
        c = 0
        for vi in @view x[lo:hi]
            if isfinite(vi)
                s += vi
                c += 1
            end
        end
        c > 0 && (out[i] = s / c)
    end
    out
end

function upsample_lidar_step(lidar_dt, mdv, vndt)
    t_src = Float64.(Dates.datetime2epochms.(lidar_dt))
    t_q = Float64.(Dates.datetime2epochms.(vndt))
    y = Float64.(mdv)
    out = fill(NaN, length(t_q))
    for i in eachindex(t_q)
        j = clamp(searchsortedlast(t_src, t_q[i]), 1, length(t_src))
        out[i] = y[j]
    end
    out
end

function one_hz_average_at_lidar_times(vndt20, x20, lidar_dt; half_window=0.5)
    out = fill(NaN, length(lidar_dt))
    t20 = Float64.(Dates.datetime2epochms.(vndt20))
    for i in eachindex(lidar_dt)
        tq = Float64(Dates.datetime2epochms(lidar_dt[i]))
        ii = findall(abs.(t20 .- tq) .<= 1000 * half_window)
        vals = filter(isfinite, x20[ii])
        isempty(vals) || (out[i] = mean(vals))
    end
    out
end

function refine_offset_20hz(win, sync_1s, Vn; ma_points=20, search_seconds=1.0)
    pad_s = ceil(Int, search_seconds) + 2
    ind20, vndt20, vn2_20 = vn_subset_20hz(win.stare_dt, Vn; pad=Second(pad_s))
    pitch_20 = Float64.(Vn[:Pitch][ind20])
    roll_20 = Float64.(Vn[:Roll][ind20])
    dt20 = dt_seconds(vndt20)
    if !isfinite(dt20) || dt20 <= 0
        dt20 = 0.05
    end

    fallback = (; final_offset_20hz=NaN, native_step=dt20, final_offset_native=sync_1s.final_offset, final_offset_round_s=round(Int, sync_1s.final_offset), vn2_1s_aligned=fill(NaN, length(win.stare_dt)), pitch_1s_aligned=fill(NaN, length(win.stare_dt)), roll_1s_aligned=fill(NaN, length(win.stare_dt)), mdv_residual_1s=fill(NaN, length(win.stare_dt)), xcorr20=nothing, vndt20, vn2_20)
    length(vndt20) < 40 && return fallback

    vn2_ma = nan_safe_moving_average(vn2_20, ma_points)
    mdv20 = upsample_lidar_step(win.stare_dt, win.mdv, vndt20)

    xcorr20 = fft_xcorr_lag(mdv20, vn2_ma; dt=dt20, center_seconds=sync_1s.final_offset, maxlag_seconds=search_seconds)
    !isfinite(xcorr20.lag_seconds) && return fallback

    final_native = round(xcorr20.lag_seconds / dt20) * dt20
    final_round_s = round(Int, final_native)

    vn2_shifted = shift_signal_linear(vn2_20, final_native; dt=dt20)
    pitch_shifted = shift_signal_linear(pitch_20, final_native; dt=dt20)
    roll_shifted = shift_signal_linear(roll_20, final_native; dt=dt20)
    vn2_1s_aligned = one_hz_average_at_lidar_times(vndt20, vn2_shifted, win.stare_dt)
    pitch_1s_aligned = one_hz_average_at_lidar_times(vndt20, pitch_shifted, win.stare_dt)
    roll_1s_aligned = one_hz_average_at_lidar_times(vndt20, roll_shifted, win.stare_dt)
    mdv_residual_1s = win.mdv .- vn2_1s_aligned

    return (; final_offset_20hz=xcorr20.lag_seconds, native_step=dt20, final_offset_native=final_native, final_offset_round_s=final_round_s, vn2_1s_aligned, pitch_1s_aligned, roll_1s_aligned, mdv_residual_1s, xcorr20, vndt20, vn2_20)
end

# UV is passed through for shared reader compatibility and later motion correction.
# The timing offsets computed here are fit only from mdv and vn2.
function process_sync_data(beams, env, Vn, UV, ic_list; ntop=80, jump_threshold_seconds=0.8, min_fine_peak=0.08, backward_windows=2, backward_tol=0.35)
    seq_results = run_sequential_offsets(beams, env, Vn, UV, ic_list; ntop=ntop, jump_threshold_seconds=jump_threshold_seconds, min_fine_peak=min_fine_peak, backward_windows=backward_windows, backward_tol=backward_tol)

    state_ref = init_stream_state()
    seq_ref20 = NamedTuple[]
    for r in seq_results
        w = extract_sync_window(beams, env, state_ref, Vn, UV, r.ic; ntop=ntop)
        s1 = coarse_and_fine_lag(w.mdv, w.vn2; prior_seconds=r.prior_offset)
        rf = refine_offset_20hz(w, s1, Vn)
        push!(seq_ref20, (; r..., lidar_dt=w.stare_dt, final_offset_20hz=rf.final_offset_20hz, final_offset_native=rf.final_offset_native, final_offset_round_s=rf.final_offset_round_s, vn2_1s_aligned=rf.vn2_1s_aligned, pitch_1s_aligned=rf.pitch_1s_aligned, roll_1s_aligned=rf.roll_1s_aligned, mdv_residual_1s=rf.mdv_residual_1s))
    end

    seq_ic = [r.ic for r in seq_ref20]
    seq_prior = [r.prior_offset for r in seq_ref20]
    seq_coarse = [r.coarse_offset for r in seq_ref20]
    seq_final = [r.final_offset for r in seq_ref20]
    seq_final_native = [r.final_offset_native for r in seq_ref20]
    seq_fine = [r.fine_residual for r in seq_ref20]

    seq_final_native_ms = round.(Int64, seq_final_native .* 1000)
    seq_jump_20hz_ms = zeros(Int64, length(seq_final_native_ms))
    if length(seq_final_native_ms) > 1
        seq_jump_20hz_ms[2:end] .= diff(seq_final_native_ms)
    end

    seq_jump_candidate = [r.jump_candidate for r in seq_ref20]
    seq_jump_accepted = [r.jump_accepted for r in seq_ref20]

    chunk_lens = [length(r.vn2_1s_aligned) for r in seq_ref20]
    chunk_record_ends = cumsum(chunk_lens)
    chunk_record_starts = vcat(1, chunk_record_ends[1:end-1] .+ 1)
    nrec = chunk_record_ends[end]

    # Schema is uniform across test/all-data runs; only ic_list size changes the output cardinality.
    return (; seq_results, seq_ref20, seq_ic, seq_prior, seq_coarse, seq_final, seq_final_native, seq_final_native_ms, seq_jump_20hz_ms, seq_fine, seq_jump_candidate, seq_jump_accepted, chunk_lens, chunk_record_starts, chunk_record_ends, nrec)
end

function save_sync_netcdf(result; nc_out=joinpath("epsilon_data", "vn_sync_offsets.nc"), time_ref_dt=DateTime(2024, 4, 29, 0, 0, 0))
    seq_ref20 = result.seq_ref20
    chunk_lens = result.chunk_lens
    chunk_record_starts = result.chunk_record_starts
    chunk_record_ends = result.chunk_record_ends
    nrec = result.nrec

    time_ref_ms = Dates.datetime2epochms(time_ref_dt)

    record_time_ms = Vector{Int64}(undef, nrec)
    record_time_dt_str = Vector{String}(undef, nrec)
    vn1s_aligned_record = Vector{Float32}(undef, nrec)
    pitch_record = Vector{Float32}(undef, nrec)
    roll_record = Vector{Float32}(undef, nrec)
    mdv_residual_record = Vector{Float32}(undef, nrec)
    record_chunk_index = Vector{Int32}(undef, nrec)

    for (i, r) in enumerate(seq_ref20)
        irec_start = chunk_record_starts[i]
        n = length(r.vn2_1s_aligned)
        for j in 1:n
            irec = irec_start + j - 1
            current_dt = r.lidar_dt[j]
            current_ms = Dates.datetime2epochms(current_dt)
            record_time_ms[irec] = Int64(current_ms - time_ref_ms)
            record_time_dt_str[irec] = string(current_dt)
            vn1s_aligned_record[irec] = Float32(r.vn2_1s_aligned[j])
            pitch_record[irec] = Float32(r.pitch_1s_aligned[j])
            roll_record[irec] = Float32(r.roll_1s_aligned[j])
            mdv_residual_record[irec] = Float32(r.mdv_residual_1s[j])
            record_chunk_index[irec] = Int32(i - 1)
        end
    end

    ds = NCDataset(nc_out, "c")
    defDim(ds, "chunk", length(seq_ref20))
    defDim(ds, "record", nrec)

    v_ic = defVar(ds, "chunk_index", Int32, ("chunk",))
    v_nrec = defVar(ds, "chunk_n_records", Int32, ("chunk",))
    v_rec_start = defVar(ds, "chunk_record_start_idx", Int32, ("chunk",))
    v_rec_end = defVar(ds, "chunk_record_end_idx", Int32, ("chunk",))
    v_t0 = defVar(ds, "chunk_start_epoch_ms", Int64, ("chunk",))
    v_t_end = defVar(ds, "chunk_end_epoch_ms", Int64, ("chunk",))
    v_off20_ms = defVar(ds, "offset_20hz_ms", Int64, ("chunk",))
    v_jump20_ms = defVar(ds, "jump_20hz_ms", Int64, ("chunk",))

    v_time = defVar(ds, "time", Int64, ("record",))
    v_time_str = defVar(ds, "time_dt", String, ("record",))
    v_chunk_idx_rec = defVar(ds, "record_chunk_index", Int32, ("record",))
    v_vn1 = defVar(ds, "vn2_1s_aligned", Float32, ("record",))
    v_pitch = defVar(ds, "pitch_degrees", Float32, ("record",))
    v_roll = defVar(ds, "roll_degrees", Float32, ("record",))
    v_res = defVar(ds, "mdv_minus_vn2_1s", Float32, ("record",))

    v_ic[:] = Int32.(result.seq_ic)
    v_nrec[:] = Int32.(chunk_lens)
    v_rec_start[:] = Int32.(chunk_record_starts)
    v_rec_end[:] = Int32.(chunk_record_ends)
    v_t0[:] = Int64.(Dates.datetime2epochms.([r.start_dt for r in seq_ref20]))
    v_t_end[:] = Int64.(Dates.datetime2epochms.([r.end_dt for r in seq_ref20]))
    v_off20_ms[:] = result.seq_final_native_ms
    v_jump20_ms[:] = result.seq_jump_20hz_ms

    v_time[:] = record_time_ms
    v_time_str[:] = record_time_dt_str
    v_chunk_idx_rec[:] = record_chunk_index
    v_vn1[:] = vn1s_aligned_record
    v_pitch[:] = pitch_record
    v_roll[:] = roll_record
    v_res[:] = mdv_residual_record

    v_time.attrib["units"] = "milliseconds since $(Dates.format(time_ref_dt, dateformat"yyyy-mm-ddTHH:MM:SS")) UTC"
    v_time.attrib["calendar"] = "standard"
    v_time.attrib["standard_name"] = "time"
    v_time.attrib["axis"] = "T"

    v_t0.attrib["long_name"] = "chunk start time in epoch milliseconds since 1970-01-01"
    v_t_end.attrib["long_name"] = "chunk end time in epoch milliseconds since 1970-01-01"
    v_off20_ms.attrib["units"] = "milliseconds"
    v_off20_ms.attrib["long_name"] = "final synchronization offset from 20 Hz native refinement"
    v_jump20_ms.attrib["units"] = "milliseconds"
    v_jump20_ms.attrib["long_name"] = "chunk-to-chunk jump in 20 Hz native offset (current minus previous)"
    v_vn1.attrib["units"] = "m s-1"
    v_vn1.attrib["long_name"] = "VectorNav VelNED2 aligned to lidar chunk times"
    v_pitch.attrib["units"] = "degrees"
    v_pitch.attrib["long_name"] = "VectorNav pitch aligned to lidar chunk times"
    v_roll.attrib["units"] = "degrees"
    v_roll.attrib["long_name"] = "VectorNav roll aligned to lidar chunk times"
    v_res.attrib["units"] = "m s-1"
    v_res.attrib["long_name"] = "mean Doppler velocity minus aligned VectorNav VelNED2"

    ds.attrib["description"] = "Lidar-VN timing offsets with 20 Hz refinement and 1-second aligned VN series (record dimension)"
    ds.attrib["Conventions"] = "CF-1.6"
    ds.attrib["time_reference_epoch"] = "$(Dates.format(time_ref_dt, dateformat"yyyy-mm-ddTHH:MM:SS")) UTC"
    close(ds)

    return nc_out
end

function read_synced_motion(nc_path=joinpath("epsilon_data", "vn_sync_offsets.nc"))
    ds = NCDataset(nc_path, "r")
    motion = (; 
        time = ds["time"][:],
        time_dt = String.(ds["time_dt"][:]),
        record_chunk_index = ds["record_chunk_index"][:],
        vn2_1s_aligned = Float64.(ds["vn2_1s_aligned"][:]),
        pitch_degrees = Float64.(ds["pitch_degrees"][:]),
        roll_degrees = Float64.(ds["roll_degrees"][:]),
        chunk_record_start_idx = ds["chunk_record_start_idx"][:],
        chunk_record_end_idx = ds["chunk_record_end_idx"][:],
    )
    close(ds)
    return motion
end

tilt_factor(pitch_degrees, roll_degrees) = cosd.(pitch_degrees) .* cosd.(roll_degrees)

function motion_correct_stare_velocity(dopplervel::AbstractVector, vn2_1s_aligned, pitch_degrees, roll_degrees)
    scale = tilt_factor(pitch_degrees, roll_degrees)
    return Float64.(dopplervel) ./ scale .- Float64.(vn2_1s_aligned)
end

function motion_correct_stare_velocity(dopplervel::AbstractMatrix, vn2_1s_aligned, pitch_degrees, roll_degrees)
    scale = reshape(tilt_factor(pitch_degrees, roll_degrees), :, 1)
    heave = reshape(Float64.(vn2_1s_aligned), :, 1)
    return Float64.(dopplervel) ./ scale .- heave
end

function read_and_motion_correct_stare(dopplervel, nc_path=joinpath("epsilon_data", "vn_sync_offsets.nc"); record_inds=:)
    motion = read_synced_motion(nc_path)
    inds = record_inds isa Colon ? collect(eachindex(motion.vn2_1s_aligned)) : record_inds
    corrected = motion_correct_stare_velocity(dopplervel, motion.vn2_1s_aligned[inds], motion.pitch_degrees[inds], motion.roll_degrees[inds])
    return (; corrected, motion=(; time=motion.time[inds], time_dt=motion.time_dt[inds], record_chunk_index=motion.record_chunk_index[inds], vn2_1s_aligned=motion.vn2_1s_aligned[inds], pitch_degrees=motion.pitch_degrees[inds], roll_degrees=motion.roll_degrees[inds]))
end

function diagnostic_single_window(beams, Env, Vn, UV, ic; ntop=80)
    single_state = init_stream_state()
    win = extract_sync_window(beams, Env, single_state, Vn, UV, ic; ntop=ntop)

    diag_array("dopplervel", win.dopplervel)
    diag_array("intensity", win.intensity)
    diag_array("Ur", win.Ur)
    diag_array("Vr", win.Vr)

    sync = coarse_and_fine_lag(win.mdv, win.vn2; prior_seconds=0.0)
    ref20 = refine_offset_20hz(win, sync, Vn)

    println("clock prior applied: ", win.lidar_clock_fast_by)
    @printf("coarse offset after iterative envelope passes = %.3f s\n", sync.coarse_offset)
    @printf("fine residual after coarse shift = %.3f s\n", sync.fine.lag_seconds)
    @printf("final offset (1 Hz stage) = %.3f s\n", sync.final_offset)
    @printf("final offset (20 Hz raw) = %.3f s\n", ref20.final_offset_20hz)
    @printf("final offset (20 Hz native step) = %.3f s [step=%.3f]\n", ref20.final_offset_native, ref20.native_step)

    fig = figure(figsize=(10, 11))
    clf()

    subplot(5, 1, 1)
    pcolormesh((beams[:time][win.ist:win.ien] .% 1) .* 60, beams[:height] ./ 1e3, pd(win.beta), cmap=ColorMap("RdYlBu_r"))
    ylabel("height (km)")
    title("Chunk $(ic): backscatter")
    colorbar()

    subplot(5, 1, 2)
    plot(win.stare_dt, win.mdv, label="mdv")
    plot(win.stare_dt, win.vn2, label="VelNED2")
    ylabel("m s^-1")
    title("Raw sync signals")
    legend()
    grid(true)

    subplot(5, 1, 3)
    plot(win.stare_dt, sync.mdv_env, label="mdv envelope")
    plot(win.stare_dt, sync.vn2_env, label="VelNED2 envelope after coarse shift")
    ylabel("arb")
    title("Envelope signals used for iterative coarse alignment")
    legend()
    grid(true)

    subplot(5, 1, 4)
    plot(win.stare_dt, sync.mdv_bp, label="mdv bandpassed")
    plot(win.stare_dt, sync.vn2_bp_coarse, label="VelNED2 after coarse shift")
    plot(win.stare_dt, sync.vn2_bp_final, label="VelNED2 after final shift", alpha=0.8)
    ylabel("arb")
    title("Real bandpassed signals used for fine FFT cross-correlation")
    legend()
    grid(true)

    subplot(5, 1, 5)
    plot(win.stare_dt, win.mdv, label="mdv")
    plot(win.stare_dt, ref20.vn2_1s_aligned, label="VN 1 s aligned")
    plot(win.stare_dt, ref20.mdv_residual_1s, label="mdv - VN", alpha=0.8)
    ylabel("m s^-1")
    title("1-second aligned VN product at lidar timing")
    legend()
    grid(true)

    tight_layout()
    return (; fig, win, sync, ref20)
end

function diagnostic_offsets(result)
    fig = figure(figsize=(5.5, 3.5))
    clf()

    ax1 = gca()
    # ax1.plot(result.seq_ic, result.seq_prior, marker="o", label="prior offset")
    # ax1.plot(result.seq_ic, result.seq_coarse, marker="o", label="coarse offset")
    ax1.plot(result.seq_ic, result.seq_final_native, marker="o", label="final offset (20 Hz native)")
    # ax1.plot(result.seq_ic, result.seq_fine, marker=".", alpha=0.5, label="fine residual")
    ax1.set_xlabel("chunk index")
    ax1.set_ylabel("seconds")
    ax1.set_title("Sequential offsets with 20 Hz native refinement")
    ax1.grid(true)

    # ax2 = ax1.twinx()
    # ax2.step(result.seq_ic, result.seq_jump_20hz_ms, where="mid", color="k", linewidth=1.4, label="jump (20 Hz native, ms)")
    # ax2.set_ylabel("jump (ms)")

    h1, l1 = ax1.get_legend_handles_labels()
    # h2, l2 = ax2.get_legend_handles_labels()
    # ax1.legend(vcat(h1, h2), vcat(l1, l2), loc="best")
    ax1.legend(h1, l1, loc="best")

    tight_layout()
    return fig
end

end