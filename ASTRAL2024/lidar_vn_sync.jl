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
using PythonPlot

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
export setup_sync_context, setup_sync_context_nc, process_sync_data, save_sync_netcdf
export write_mdv_sync_pass2!, write_daily_mdv_vn2!
export diagnostic_single_window, diagnostic_offsets, diagnostic_example_chunks
export read_synced_motion, motion_correct_stare_velocity, read_and_motion_correct_stare

const NOISE_THR = 1.03
const BANDPASS_PERIOD = (5.0, 20.0)
const BANDPASS_HZ = (1 / BANDPASS_PERIOD[2], 1 / BANDPASS_PERIOD[1])
const TIMESTEP = 1.02
const RANGEGATE = 24.0
const OFFSET_SENTINEL_S = -9999.0
const OFFSET_SENTINEL_MS = Int64(-9_999_000)
const NAN_OFFSET_LOG_PATH = joinpath("epsilon_data", "nan_offset_chunks.log")

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
    ncdir = joinpath(lidarstemdir, "netcdf_stare")
    starefiles = sort(filter(f -> startswith(f, "Stare_116_") && endswith(f, ".nc"), readdir(ncdir)))
    ff = joinpath.(ncdir, starefiles)
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

"""
ensure_chunk_loaded_nc!(beams, env, state, ic; nc_dir)

Drop-in replacement for `ensure_chunk_loaded!` that reads pre-converted
NetCDF files from `nc_dir` (default `data/netcdf_stare/`) instead of
parsing raw HPL text files.  All other sync machinery is unchanged.
"""
function ensure_chunk_loaded_nc!(beams, env, state, ic;
                                  nc_dir::AbstractString = joinpath("data", "netcdf_stare"))
    ist = env.ists[ic]
    ien = env.iens[ic]
    ifile_needed = findlast(env.bigind_file_starts .<= ist)
    isnothing(ifile_needed) && error("No file contains ist=$(ist)")

    if state[:ifile_loaded] == 0
        state[:ifile_loaded] = ifile_needed - 1
        state[:bigind_file_end] = 0
    end

    nz = size(beams[:height][:], 1)

    while (state[:ifile_loaded] < ifile_needed) || (ien > state[:bigind_file_end])
        state[:ifile_loaded] += 1
        state[:ifile_loaded] > length(env.ff) &&
            error("Ran out of lidar files while loading chunk ic=$(ic)")

        ifile = state[:ifile_loaded]
        bb    = env.bigind_file_starts[ifile]:env.bigind_file_ends[ifile]

        # env.ff already contains the NC path directly
        nc_path = env.ff[ifile]

        NCDatasets.NCDataset(nc_path, "r") do ds
            t_nc = ds["time"][:]
            n = length(t_nc)
            n != length(bb) && error("NetCDF beam count mismatch for $(basename(nc_path)): nc=$(n), expected=$(length(bb))")

            # CF-decoded time may come back as DateTime; convert to HPL-style decimal hours.
            t_hpl_hours = Vector{Float32}(undef, n)
            for i in eachindex(t_nc)
                ti = t_nc[i]
                if ti isa Dates.TimeType
                    t_hpl_hours[i] = Float32(
                        Dates.hour(ti) +
                        Dates.minute(ti) / 60 +
                        Dates.second(ti) / 3600 +
                        Dates.millisecond(ti) / 3_600_000
                    )
                elseif ismissing(ti)
                    t_hpl_hours[i] = Float32(NaN)
                else
                    t_hpl_hours[i] = Float32(ti)
                end
            end

            ng = min(size(ds["dopplervel"], 2), nz)

            setindex!(beams[:time],      t_hpl_hours,                    bb)
            setindex!(beams[:azimuth],   Float32.(ds["azimuth"][1:n]),   bb)
            setindex!(beams[:elevangle], Float32.(ds["elevangle"][1:n]), bb)
            setindex!(beams[:pitch],     Float32.(ds["pitch"][1:n]),     bb)
            setindex!(beams[:roll],      Float32.(ds["roll"][1:n]),      bb)
            beams[:height][1:ng] .= Float32.(ds["height"][1:ng])
            setindex!(beams[:dopplervel], Float32.(ds["dopplervel"][1:n, 1:ng]), bb, 1:ng)
            setindex!(beams[:intensity],  Float32.(ds["intensity"][1:n,  1:ng]), bb, 1:ng)
            setindex!(beams[:beta],       Float32.(ds["beta"][1:n,       1:ng]), bb, 1:ng)
        end

        state[:bigind_file_end] = env.bigind_file_ends[ifile]
    end

    return (; ifile=ifile_needed, ist, ien)
end

"""
setup_sync_context_nc(; nc_dir, ...)

Like `setup_sync_context` but patches `extract_sync_window` to load beam
data from pre-converted NetCDF files rather than raw HPL text.

The returned context has an extra field `nc_dir` and a `load_nc!` closure
so callers can swap in the NC loader by passing `load_fn=ctx.load_nc!` to
`ensure_chunk_loaded_nc!`.  In practice the notebook just uses the
`beams_nc` variant of `extract_sync_window` defined below.
"""
function setup_sync_context_nc(;
        lidarstemdir = "./data",
        uv_path      = joinpath("data/netcdf", "ekamsat_lidar_uv_20240428-20240613.nc"),
        nc_dir       = joinpath("data", "netcdf_stare"),
        nx           = 4000,
        nz           = 80)

    Env = load_lidar_indices_and_files(; lidarstemdir=lidarstemdir)
    Vn  = read_vecnav_dict()
    UV  = NCDatasets.NCDataset(uv_path)

    dtime_st = Env.dtime[Env.ists]
    dtime_en = Env.dtime[Env.iens]
    icvn     = findfirst(dtime_st .>= Vn[:vndt][1]):findlast(dtime_en .<= Vn[:vndt][end])
    beams    = init_periodic_beams(nx, nz)

    return (; Env, Vn, UV, beams, icvn, nz, nx, nc_dir)
end

function chunk_lidar_datetimes(dt_chunk, beams, ist, ien)
    base_date = Date(dt_chunk)
    stare_dt_raw = @. DateTime(base_date) + Millisecond(round(Int64, beams[:time][ist:ien] * 3_600_000))
    # Detect midnight crossing: if any timestamp steps backward, advance
    # all subsequent timestamps by one day (one crossing per chunk is enough).
    for i in 2:length(stare_dt_raw)
        if stare_dt_raw[i] < stare_dt_raw[i-1]
            stare_dt_raw[i:end] .+= Day(1)
            break
        end
    end
    lidar_clock_fast_by = Millisecond(round(Int64, 1_000 * fit_offset(stare_dt_raw[1])))
    stare_dt = stare_dt_raw .- lidar_clock_fast_by
    return stare_dt_raw, stare_dt, lidar_clock_fast_by
end

"""
extract_sync_window
Subsets Vn, and lidar beams.
Calls ensure_chunk_loaded to update `beams` in-place.
UV, not used here, is passed to read_stare_chunk. 
The sync offset will depend only on mdv and vn2, not on Ur,Vr.
"""
function extract_sync_window(beams, env, state, Vn, UV, ic; ntop=80, nc_dir="data/netcdf_stare")
    if isnothing(nc_dir)
        info = ensure_chunk_loaded!(beams, env, state, ic)
    else
        info = ensure_chunk_loaded_nc!(beams, env, state, ic; nc_dir=nc_dir)
    end
    dt_arg = info.ist <= length(env.dtime) ? env.dtime[info.ist] : env.dtime[end]
    dopplervel, pitch, roll, vn0, vn1, vn2, Ur, Vr, mdv_builtin = read_stare_chunk(dt_arg, beams, Vn, UV, info.ist, info.ien, ntop)
    stare_dt_raw, stare_dt, lidar_clock_fast_by = chunk_lidar_datetimes(env.dtime[info.ist], beams, info.ist, info.ien)
    intensity = beams[:intensity][info.ist:info.ien, 1:ntop]
    beta = beams[:beta][info.ist:info.ien, 1:ntop]
    mdv, goodlevels = valid_chunk_mdv(intensity, dopplervel)
    _, vndt20_chunk, vn2_20_chunk = vn_subset_20hz(stare_dt, Vn; pad=Second(2))
    vndt20_chunk, vn2_20_chunk, _ = remove_vn_in_lidar_gaps(stare_dt, vndt20_chunk, vn2_20_chunk)
    vn2_xcorr = one_hz_average_at_lidar_times(vndt20_chunk, vn2_20_chunk, stare_dt)
    vn_coverage = vn_coverage_fraction(stare_dt, Vn)
    vn2_xcorr_nan_frac = count(!isfinite, vn2_xcorr) / max(length(vn2_xcorr), 1)
    return (; info..., stare_dt_raw, stare_dt, lidar_clock_fast_by, dopplervel, intensity, beta, pitch, roll, vn0, vn1, vn2, vn2_xcorr, vndt20_chunk, Ur, Vr, mdv, mdv_builtin, goodlevels, vn_coverage, vn2_xcorr_nan_frac)
end

"""
    remove_vn_in_lidar_gaps(stare_dt, vndt, vn2; gap_threshold_s=1.5)

Remove VN samples whose timestamps fall within lidar gap intervals (periods where
the lidar has no beams, detected as consecutive beam times > gap_threshold_s apart,
e.g. at file rollovers). This prevents `one_hz_average_at_lidar_times` from mixing
VN data across the gap when computing per-beam averages near the gap boundary.

Returns `(vndt_filtered, vn2_filtered, keep)` where `keep` is a BitVector.
"""
function remove_vn_in_lidar_gaps(stare_dt, vndt, vn2; gap_threshold_s=1.5)
    keep = trues(length(vndt))
    length(stare_dt) < 2 && return vndt, vn2, keep
    beam_diffs_ms = Float64.(Dates.value.(Millisecond.(diff(stare_dt))))
    gap_mask = beam_diffs_ms .> (gap_threshold_s * 1000)
    any(gap_mask) || return vndt, vn2, keep
    t_vn_ms = Float64.(Dates.datetime2epochms.(vndt))
    for ig in findall(gap_mask)
        t_gap_start = Float64(Dates.datetime2epochms(stare_dt[ig]))
        t_gap_end   = Float64(Dates.datetime2epochms(stare_dt[ig + 1]))
        for j in eachindex(t_vn_ms)
            if t_vn_ms[j] > t_gap_start && t_vn_ms[j] < t_gap_end
                keep[j] = false
            end
        end
    end
    return vndt[keep], vn2[keep], keep
end

function finite_overlap(x, y)
    n = min(length(x), length(y))
    xv = x[1:n]
    yv = y[1:n]
    keep = isfinite.(xv) .& isfinite.(yv)
    return xv[keep], yv[keep], keep
end

function fill_short_nan_gaps(x; max_gap=3)
    y = Float64.(x)
    n = length(y)
    i = 1
    while i <= n
        if isfinite(y[i])
            i += 1
            continue
        end
        i0 = i
        while i <= n && !isfinite(y[i])
            i += 1
        end
        i1 = i - 1
        gap = i1 - i0 + 1
        if gap <= max_gap && i0 > 1 && i <= n && isfinite(y[i0 - 1]) && isfinite(y[i])
            yl = y[i0 - 1]
            yr = y[i]
            for k in 0:(gap - 1)
                w = (k + 1) / (gap + 1)
                y[i0 + k] = (1 - w) * yl + w * yr
            end
        end
    end
    y
end

"""
    vn_coverage_fraction(stare_dt, Vn)

Return the fraction of the lidar chunk's time span that is covered by VectorNav 20 Hz data.
A value < 0.5 means fewer than half the expected samples are present; the chunk is treated as
having insufficient VN coverage and no timing offset is computed.
"""
function vn_coverage_fraction(stare_dt, Vn)
    t0 = stare_dt[1]
    t1 = stare_dt[end]
    span_s = Dates.value(Millisecond(t1 - t0)) / 1000.0
    span_s <= 0.0 && return 0.0
    n_vn = count(t0 .<= Vn[:vndt] .<= t1)
    expected = span_s * 20.0   # nominal 20 Hz
    return n_vn / max(expected, 1.0)
end

offset_or_sentinel(x; sentinel=-9999.0) = isfinite(x) ? x : sentinel
offset_to_ms(x; sentinel_ms=Int64(-9_999_000)) = isfinite(x) ? round(Int64, x * 1000) : sentinel_ms

function append_nan_offset_log(log_path; ic, ist, ien, start_dt, end_dt, reason, prior_offset, coarse_offset, fine_residual, final_offset_1hz, final_offset_20hz, final_offset_native, nvalid_mdv, nvalid_vn2)
    mkpath(dirname(log_path))
    write_header = !isfile(log_path)
    open(log_path, "a") do io
        if write_header
            println(io, "run_utc,ic,ist,ien,start_dt,end_dt,reason,prior_offset_s,coarse_offset_s,fine_residual_s,final_offset_1hz_s,final_offset_20hz_s,final_offset_native_s,nvalid_mdv,nvalid_vn2")
        end
        println(io, string(Dates.now(Dates.UTC), ",", ic, ",", ist, ",", ien, ",", start_dt, ",", end_dt, ",", reason, ",", prior_offset, ",", coarse_offset, ",", fine_residual, ",", final_offset_1hz, ",", final_offset_20hz, ",", final_offset_native, ",", nvalid_mdv, ",", nvalid_vn2))
    end
    nothing
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

function coarse_and_fine_lag(mdv, vn2; dt=TIMESTEP, prior_seconds=0.0, coarse_search_seconds=12.0, max_passes=4, fine_search_seconds=4.0, max_gap_samples=3, max_vn_nan_frac=0.15)
    vn_nan_frac = count(!isfinite, Float64.(vn2)) / max(length(vn2), 1)
    if vn_nan_frac > max_vn_nan_frac
        nan_arr = fill(NaN, length(vn2))
        nan_fine = (lag_seconds=NaN, peak=NaN, peak_norm=NaN, lags=Float64[], corr=Float64[])
        return (; prior_seconds, mdv_bp=nan_arr, vn2_bp=nan_arr, mdv_env=nan_arr, vn2_env=nan_arr, coarse_history=NamedTuple[], coarse_offset=NaN, fine=nan_fine, final_offset=NaN, vn2_bp_coarse=nan_arr, vn2_bp_final=nan_arr, vn2_env_final=nan_arr)
    end
    mdv_clean = fill_short_nan_gaps(mdv; max_gap=max_gap_samples)
    vn2_clean = fill_short_nan_gaps(vn2; max_gap=max_gap_samples)
    mdv_bp = fft_bandpass(mdv_clean; dt=dt)
    vn2_bp = fft_bandpass(vn2_clean; dt=dt)

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

function backward_jump_robustness(beams, env, state, Vn, UV, history, post_offset; ntop=80, n_back=2, max_gap_samples=3, nc_dir=nothing)
    isempty(history) && return NaN
    ntest = min(n_back, length(history))
    deltas = Float64[]
    for k in 1:ntest
        h = history[end - k + 1]
        win_prev = extract_sync_window(beams, env, state, Vn, UV, h.ic; ntop=ntop, nc_dir=nc_dir)
        back_sync = coarse_and_fine_lag(win_prev.mdv, win_prev.vn2_xcorr; prior_seconds=post_offset, max_gap_samples=max_gap_samples)
        if isfinite(back_sync.final_offset) && isfinite(h.final_offset)
            push!(deltas, abs(back_sync.final_offset - h.final_offset))
        end
    end
    isempty(deltas) ? NaN : median(deltas)
end

function run_sequential_offsets(beams, env, Vn, UV, ic_list; ntop=80, jump_threshold_seconds=0.8, min_fine_peak=0.08, backward_windows=2, backward_tol=0.35, max_gap_samples=3, min_accept_peak_norm=0.08, fallback_search_seconds=60.0, wide_search_seconds=30.0, reset_prior_on_file_boundary=false, reset_prior_on_rejected_sync=false, nc_dir="data/netcdf_stare")
    state = init_stream_state()
    history = NamedTuple[]
    last_valid_offset = NaN
    prev_failed = false
    for ic in ic_list
        win = extract_sync_window(beams, env, state, Vn, UV, ic; ntop=ntop, nc_dir=nc_dir)
        start_dt = win.stare_dt[1]
        prior = prior_from_history(history, start_dt)

        prior_reset = "none"
        # If previous chunk failed, use last valid offset as prior for this chunk
        if prev_failed && isfinite(last_valid_offset)
            prior = (; prior..., prior_seconds=last_valid_offset)
            prior_reset = "carry_forward"
        elseif !isempty(history)
            if reset_prior_on_file_boundary && hasproperty(history[end], :ifile) && history[end].ifile != win.ifile
                prior = (; prior..., prior_seconds=0.0)
                prior_reset = "file_boundary"
            elseif reset_prior_on_rejected_sync && !get(history[end], :accepted_sync, true)
                prior = (; prior..., prior_seconds=0.0)
                prior_reset = "rejected_sync"
            end
        end

        if win.vn_coverage < 0.5 || win.vn2_xcorr_nan_frac > 0.15
            prev_failed = true
            push!(history, (; ic, ifile=win.ifile, start_dt=win.stare_dt[1], end_dt=win.stare_dt[end], parent_start=floor(win.stare_dt[1], Hour), prior_offset=prior.prior_seconds, previous_offset=prior.previous_offset, recent_offset=prior.recent_offset, parent_offset=prior.parent_offset, prior_reset, coarse_offset=NaN, fine_residual=NaN, final_offset=NaN, coarse_peak=NaN, fine_peak=NaN, goodlevels=length(win.goodlevels), jump_candidate=false, jump_robust=false, jump_accepted=false, jump_backward_metric=NaN, fallback_used=false, accepted_sync=false))
            continue
        end
        sync = coarse_and_fine_lag(win.mdv, win.vn2_xcorr; prior_seconds=prior.prior_seconds, max_gap_samples=max_gap_samples)

        # Always compare against a prior-independent wide search (±wide_search_seconds from 0).
        # This prevents a bad prior from trapping the iterative coarse stage at a local-max
        # envelope peak. The prior-guided result is kept only if it has better fine_peak_norm.
        fallback_used = false
        sync_wide = coarse_and_fine_lag(win.mdv, win.vn2_xcorr;
            prior_seconds=0.0, coarse_search_seconds=wide_search_seconds,
            max_gap_samples=max_gap_samples)
        if isfinite(sync_wide.fine.peak_norm) &&
           (!isfinite(sync.fine.peak_norm) || sync_wide.fine.peak_norm > sync.fine.peak_norm)
            sync = sync_wide
            fallback_used = true
        end

        jump_candidate = is_jump_candidate(prior.prior_seconds, sync.final_offset, sync.fine.peak_norm; jump_threshold_seconds=jump_threshold_seconds, min_fine_peak=min_fine_peak)

        jump_backward_metric = NaN
        jump_robust = false
        jump_accepted = false

        if jump_candidate && !isempty(history)
            post_sync = coarse_and_fine_lag(win.mdv, win.vn2_xcorr; prior_seconds=sync.final_offset, max_gap_samples=max_gap_samples)
            post_offset = post_sync.final_offset
            jump_backward_metric = backward_jump_robustness(beams, env, state, Vn, UV, history, post_offset; ntop=ntop, n_back=backward_windows, max_gap_samples=max_gap_samples, nc_dir=nc_dir)
            jump_robust = isfinite(jump_backward_metric) && jump_backward_metric <= backward_tol
            if jump_robust
                sync = post_sync
                jump_accepted = true
            else
                sync = coarse_and_fine_lag(win.mdv, win.vn2_xcorr; prior_seconds=prior.prior_seconds, fine_search_seconds=2.0, max_gap_samples=max_gap_samples)
            end
        end

        accepted_sync = isfinite(sync.fine.peak_norm) && sync.fine.peak_norm >= min_accept_peak_norm
        stored_coarse_offset = accepted_sync ? sync.coarse_offset : NaN
        stored_fine_residual = accepted_sync ? sync.fine.lag_seconds : NaN
        stored_final_offset = accepted_sync ? sync.final_offset : NaN

        if accepted_sync && isfinite(stored_final_offset)
            last_valid_offset = stored_final_offset
            prev_failed = false
        else
            prev_failed = true
        end

        push!(history, (; ic, ifile=win.ifile, start_dt=win.stare_dt[1], end_dt=win.stare_dt[end], parent_start=floor(win.stare_dt[1], Hour), prior_offset=prior.prior_seconds, previous_offset=prior.previous_offset, recent_offset=prior.recent_offset, parent_offset=prior.parent_offset, prior_reset, coarse_offset=stored_coarse_offset, fine_residual=stored_fine_residual, final_offset=stored_final_offset, coarse_peak=isempty(sync.coarse_history) ? NaN : sync.coarse_history[end].peak_norm, fine_peak=sync.fine.peak_norm, goodlevels=length(win.goodlevels), jump_candidate, jump_robust, jump_accepted, jump_backward_metric, fallback_used, accepted_sync))
    end
    history
end

function vn_subset_20hz(stare_dt, Vn; pad=Second(2), offset_seconds=0.0)
    dt_off = Millisecond(round(Int64, 1_000 * offset_seconds))
    t0 = stare_dt[1] + dt_off
    t1 = stare_dt[end] + dt_off
    ind = findall(t0 - pad .<= Vn[:vndt] .<= t1 + pad)
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

function one_hz_average_at_lidar_times(vndt20, x20, lidar_dt; half_window=0.5, query_offset_seconds=0.0)
    out = fill(NaN, length(lidar_dt))
    t20 = Float64.(Dates.datetime2epochms.(vndt20))
    qoff_ms = 1_000 * query_offset_seconds
    for i in eachindex(lidar_dt)
        tq = Float64(Dates.datetime2epochms(lidar_dt[i])) + qoff_ms
        ii = findall(abs.(t20 .- tq) .<= 1000 * half_window)
        vals = filter(isfinite, x20[ii])
        isempty(vals) || (out[i] = mean(vals))
    end
    out
end

function refine_offset_20hz(win, sync_1s, Vn; ma_points=5, search_seconds=0.6, max_gap_samples=3)
    prior_offset = isfinite(sync_1s.prior_seconds) ? sync_1s.prior_seconds : 0.0
    offset_base_1hz = isfinite(sync_1s.final_offset) ? sync_1s.final_offset : prior_offset
    pad_s = ceil(Int, search_seconds) + 2
    ind20, vndt20, vn2_20 = vn_subset_20hz(win.stare_dt, Vn; pad=Second(pad_s), offset_seconds=-offset_base_1hz)
    vndt20, vn2_20, vn_keep = remove_vn_in_lidar_gaps(win.stare_dt, vndt20, vn2_20)
    ind20 = ind20[vn_keep]
    pitch_20 = Float64.(Vn[:Pitch][ind20])
    roll_20 = Float64.(Vn[:Roll][ind20])
    base_offset = isfinite(sync_1s.final_offset) ? sync_1s.final_offset : (isfinite(sync_1s.prior_seconds) ? sync_1s.prior_seconds : -9999.0)
    dt20 = dt_seconds(vndt20)
    if !isfinite(dt20) || dt20 <= 0
        dt20 = 0.05
    end

    # Fallback means "no additional 20 Hz correction"; keep residual at 0 s.
    fallback = (; final_offset_20hz=0.0, native_step=dt20, final_offset_native=0.0, final_offset_round_s=0, vn2_1s_aligned=fill(NaN, length(win.stare_dt)), pitch_1s_aligned=fill(NaN, length(win.stare_dt)), roll_1s_aligned=fill(NaN, length(win.stare_dt)), mdv_residual_1s=fill(NaN, length(win.stare_dt)), xcorr20=nothing, vndt20, vn2_20)
    length(vndt20) < 40 && return fallback

    dt1 = dt_seconds(win.stare_dt)
    (!isfinite(dt1) || dt1 <= 0) && (dt1 = TIMESTEP)
    
    vn2_20_ma = nan_safe_moving_average(vn2_20, ma_points)

    # Evaluate residual lags at native 20 Hz cadence (0.05 s).
    search_step = 0.05
    search_steps = collect(-search_seconds:search_step:search_seconds)
    best_corr = -Inf
    best_offset = 0.0
    
    mdv1s = fill_short_nan_gaps(win.mdv; max_gap=max_gap_samples)
    
    for offset_candidate in search_steps
        vn2_20_shifted = shift_signal_linear(vn2_20_ma, offset_candidate; dt=dt20)
        vn2_1s = one_hz_average_at_lidar_times(vndt20, vn2_20_shifted, win.stare_dt; query_offset_seconds=-offset_base_1hz)
        
        # Compute correlation at 1 Hz (finite_overlap_corr handles NaN safely; do not fill VN gaps)
        corr = finite_overlap_corr(mdv1s, vn2_1s)
        
        if isfinite(corr) && corr > best_corr
            best_corr = corr
            best_offset = offset_candidate
        end
    end
    
    # Apply best offset to all 20 Hz signals for output
    vn2_shifted = shift_signal_linear(vn2_20_ma, best_offset; dt=dt20)
    pitch_shifted = shift_signal_linear(pitch_20, best_offset; dt=dt20)
    roll_shifted = shift_signal_linear(roll_20, best_offset; dt=dt20)
    vn2_1s_aligned = one_hz_average_at_lidar_times(vndt20, vn2_shifted, win.stare_dt; query_offset_seconds=-offset_base_1hz)
    pitch_1s_aligned = one_hz_average_at_lidar_times(vndt20, pitch_shifted, win.stare_dt; query_offset_seconds=-offset_base_1hz)
    roll_1s_aligned = one_hz_average_at_lidar_times(vndt20, roll_shifted, win.stare_dt; query_offset_seconds=-offset_base_1hz)
    mdv_residual_1s = win.mdv .- vn2_1s_aligned

    final_round_s = round(Int, best_offset)
    
    # Construct minimal xcorr struct for compatibility
    xcorr20 = (; lag_seconds=best_offset, corr=best_corr)

    return (; final_offset_20hz=best_offset, native_step=dt20, final_offset_native=best_offset, final_offset_round_s=final_round_s, vn2_1s_aligned, pitch_1s_aligned, roll_1s_aligned, mdv_residual_1s, xcorr20, vndt20, vn2_20)
end

# Helper: Pearson correlation for finite overlaps
function finite_overlap_corr(x, y)
    xv, yv, _ = finite_overlap(Float64.(x), Float64.(y))
    length(xv) < 4 && return NaN
    sx = std(xv)
    sy = std(yv)
    (!isfinite(sx) || !isfinite(sy) || sx == 0 || sy == 0) && return NaN
    return cor(xv, yv)
end

# UV is passed through for shared reader compatibility and later motion correction.
# The timing offsets computed here are fit only from mdv and vn2.
function process_sync_data(beams, env, Vn, UV, ic_list; 
    ntop=80, jump_threshold_seconds=0.8, min_fine_peak=0.08, 
    backward_windows=2, backward_tol=0.35, max_gap_samples=3,
    min_accept_peak_norm=0.08, fallback_search_seconds=60.0, wide_search_seconds=30.0,
    reset_prior_on_file_boundary=false, reset_prior_on_rejected_sync=false,
    nc_dir="data/netcdf_stare",
    offset_sentinel_seconds=-9999.0, offset_sentinel_ms=Int64(-9_999_000),
    nan_log_path=joinpath("epsilon_data", "nan_offset_chunks.log"), 
    reset_nan_log=true, 
    iter_log_path=joinpath("epsilon_data", "vn_log_$(Dates.format(Dates.now(), dateformat"yyyymmdd_HHMMSS")).txt"), 
    reset_iter_log=true )
    
    if reset_nan_log && isfile(nan_log_path)
        rm(nan_log_path)
    end
    if reset_iter_log && isfile(iter_log_path)
        rm(iter_log_path)
    end
    mkpath(dirname(iter_log_path))

    seq_results = run_sequential_offsets(beams, env, Vn, UV, ic_list; ntop=ntop, jump_threshold_seconds=jump_threshold_seconds, min_fine_peak=min_fine_peak, backward_windows=backward_windows, backward_tol=backward_tol, max_gap_samples=max_gap_samples, min_accept_peak_norm=min_accept_peak_norm, fallback_search_seconds=fallback_search_seconds, wide_search_seconds=wide_search_seconds, reset_prior_on_file_boundary=reset_prior_on_file_boundary, reset_prior_on_rejected_sync=reset_prior_on_rejected_sync, nc_dir=nc_dir)

    state_ref = init_stream_state()
    seq_ref20 = NamedTuple[]
    open(iter_log_path, "a") do iter_io

        if filesize(iter_log_path) == 0
            println(iter_io, "iter,ic,timestamp_utc,status,message")
            flush(iter_io)
        end

        for (iter, r) in enumerate(seq_results)
            println(iter_io, string(iter, ",", r.ic, ",", Dates.now(Dates.UTC), ",start,"))
            flush(iter_io)
            status = "ok"
            message = ""

            try
                w = extract_sync_window(beams, env, state_ref, Vn, UV, r.ic; ntop=ntop, nc_dir=nc_dir)
                if w.vn_coverage < 0.5 || w.vn2_xcorr_nan_frac > 0.15
                    status = "sentinel"
                    message = w.vn_coverage < 0.5 ? "insufficient_vn_coverage" : "high_vn2_xcorr_nan_frac"
                    append_nan_offset_log(nan_log_path;
                        ic=r.ic,
                        ist=w.ist,
                        ien=w.ien,
                        start_dt=w.stare_dt[1],
                        end_dt=w.stare_dt[end],
                        reason=message,
                        prior_offset=r.prior_offset,
                        coarse_offset=NaN,
                        fine_residual=NaN,
                        final_offset_1hz=NaN,
                        final_offset_20hz=NaN,
                        final_offset_native=NaN,
                        nvalid_mdv=count(isfinite, w.mdv),
                        nvalid_vn2=count(isfinite, w.vn2_xcorr),
                    )
                    push!(seq_ref20, (; r..., lidar_dt=w.stare_dt, final_offset_20hz=offset_sentinel_seconds, final_offset_native=offset_sentinel_seconds, final_offset_round_s=round(Int, offset_sentinel_seconds), vn2_1s_aligned=fill(NaN, length(w.stare_dt)), pitch_1s_aligned=fill(NaN, length(w.stare_dt)), roll_1s_aligned=fill(NaN, length(w.stare_dt)), mdv_residual_1s=fill(NaN, length(w.stare_dt))))
                elseif !get(r, :accepted_sync, true)
                    status = "sentinel"
                    message = "rejected_low_confidence_sync"
                    append_nan_offset_log(nan_log_path;
                        ic=r.ic,
                        ist=w.ist,
                        ien=w.ien,
                        start_dt=w.stare_dt[1],
                        end_dt=w.stare_dt[end],
                        reason=message,
                        prior_offset=r.prior_offset,
                        coarse_offset=r.coarse_offset,
                        fine_residual=r.fine_residual,
                        final_offset_1hz=r.final_offset,
                        final_offset_20hz=NaN,
                        final_offset_native=NaN,
                        nvalid_mdv=count(isfinite, w.mdv),
                        nvalid_vn2=count(isfinite, w.vn2_xcorr),
                    )
                    # Fallback when refine_offset_20hz fails: use 1 Hz result directly
                    final_native_fallback = isfinite(r.final_offset) ? r.final_offset : offset_sentinel_seconds
                    push!(seq_ref20, (; r..., lidar_dt=w.stare_dt, final_offset_20hz=offset_sentinel_seconds, final_offset_native=final_native_fallback, final_offset_round_s=(isfinite(final_native_fallback) ? round(Int, final_native_fallback) : round(Int, offset_sentinel_seconds)), vn2_1s_aligned=fill(NaN, length(w.stare_dt)), pitch_1s_aligned=fill(NaN, length(w.stare_dt)), roll_1s_aligned=fill(NaN, length(w.stare_dt)), mdv_residual_1s=fill(NaN, length(w.stare_dt))))
                else
                    s1 = coarse_and_fine_lag(w.mdv, w.vn2_xcorr; prior_seconds=r.prior_offset, max_gap_samples=max_gap_samples)
                    # Always compare against a prior-independent wide search.
                    s1_wide = coarse_and_fine_lag(w.mdv, w.vn2_xcorr;
                        prior_seconds=0.0, coarse_search_seconds=wide_search_seconds,
                        max_gap_samples=max_gap_samples)
                    if isfinite(s1_wide.fine.peak_norm) &&
                       (!isfinite(s1.fine.peak_norm) || s1_wide.fine.peak_norm > s1.fine.peak_norm)
                        s1 = s1_wide
                    end
                    rf = refine_offset_20hz(w, s1, Vn; max_gap_samples=max_gap_samples)

                        final_20 = offset_or_sentinel(rf.final_offset_20hz; sentinel=offset_sentinel_seconds)
                        final_native_resid = offset_or_sentinel(rf.final_offset_native; sentinel=offset_sentinel_seconds)
                        # final_offset_native from refine_offset_20hz is the 20 Hz residual (no prior in it).
                        # s1.final_offset is the 1 Hz result which already includes prior.
                        # Correct composition: s1.final_offset + residual_20hz (no double-counting of prior)
                        final_native = isfinite(final_native_resid) && isfinite(s1.final_offset) ? s1.final_offset + final_native_resid : offset_sentinel_seconds
                        final_round = isfinite(final_native) && final_native != offset_sentinel_seconds ? round(Int, final_native) : round(Int, offset_sentinel_seconds)

                    if !isfinite(rf.final_offset_20hz) || !isfinite(rf.final_offset_native)
                        reason = !isfinite(rf.final_offset_native) ? "nan_native_offset" : "nan_20hz_offset"
                        append_nan_offset_log(nan_log_path;
                            ic=r.ic,
                            ist=w.ist,
                            ien=w.ien,
                            start_dt=w.stare_dt[1],
                            end_dt=w.stare_dt[end],
                            reason=reason,
                            prior_offset=s1.prior_seconds,
                            coarse_offset=s1.coarse_offset,
                            fine_residual=s1.fine.lag_seconds,
                            final_offset_1hz=s1.final_offset,
                            final_offset_20hz=rf.final_offset_20hz,
                            final_offset_native=rf.final_offset_native,
                            nvalid_mdv=count(isfinite, w.mdv),
                            nvalid_vn2=count(isfinite, w.vn2_xcorr),
                        )
                        status = "sentinel"
                        message = reason
                    end

                    push!(seq_ref20, (; r..., lidar_dt=w.stare_dt, final_offset_20hz=final_20, final_offset_native=final_native, final_offset_round_s=final_round, vn2_1s_aligned=rf.vn2_1s_aligned, pitch_1s_aligned=rf.pitch_1s_aligned, roll_1s_aligned=rf.roll_1s_aligned, mdv_residual_1s=rf.mdv_residual_1s))
                end
            catch err
                status = "error"
                message = replace(sprint(showerror, err), '\n' => " | ")
                println(iter_io, string(iter, ",", r.ic, ",", Dates.now(Dates.UTC), ",", status, ",", message))
                flush(iter_io)
                rethrow(err)
            end

            println(iter_io, string(iter, ",", r.ic, ",", Dates.now(Dates.UTC), ",", status, ",", message))
            flush(iter_io)
        end
    end

    seq_ic = [r.ic for r in seq_ref20]
    seq_prior = [r.prior_offset for r in seq_ref20]
    seq_coarse = [r.coarse_offset for r in seq_ref20]
    seq_final = [r.final_offset for r in seq_ref20]
    seq_final_native = [r.final_offset_native for r in seq_ref20]
    seq_fine = [r.fine_residual for r in seq_ref20]

    seq_final_native_ms = [offset_to_ms(x; sentinel_ms=offset_sentinel_ms) for x in seq_final_native]
    seq_jump_20hz_ms = zeros(Int64, length(seq_final_native_ms))
    if length(seq_final_native_ms) > 1
        for i in 2:length(seq_final_native_ms)
            if seq_final_native_ms[i] == offset_sentinel_ms || seq_final_native_ms[i - 1] == offset_sentinel_ms
                seq_jump_20hz_ms[i] = offset_sentinel_ms
            else
                seq_jump_20hz_ms[i] = seq_final_native_ms[i] - seq_final_native_ms[i - 1]
            end
        end
    end

    seq_jump_candidate = [r.jump_candidate for r in seq_ref20]
    seq_jump_accepted = [r.jump_accepted for r in seq_ref20]
    seq_fallback_used = [r.fallback_used for r in seq_ref20]
    seq_accepted_sync = [get(r, :accepted_sync, true) for r in seq_ref20]

    chunk_lens = [length(r.vn2_1s_aligned) for r in seq_ref20]
    chunk_record_ends = cumsum(chunk_lens)
    chunk_record_starts = vcat(1, chunk_record_ends[1:end-1] .+ 1)
    nrec = chunk_record_ends[end]

    # Schema is uniform across test/all-data runs; only ic_list size changes the output cardinality.
    return (; seq_results, seq_ref20, seq_ic, seq_prior, seq_coarse, seq_final, seq_final_native, seq_final_native_ms, seq_jump_20hz_ms, seq_fine, seq_jump_candidate, seq_jump_accepted, seq_fallback_used, seq_accepted_sync, chunk_lens, chunk_record_starts, chunk_record_ends, nrec)
end

"""
write_mdv_sync_pass2!(; nc_dir, lidarstemdir, nx, nz, ntop, ic_list, overwrite, log_every)

Second pass that appends `mdv_sync(time)` to each hourly stare NetCDF file.
`mdv_sync` is computed chunk-by-chunk using the same `extract_sync_window`
path as the sync loop, so it is exactly the same signal used by timing sync.

This function leaves existing `mdv_snr_mean` unchanged.
"""
function write_mdv_sync_pass2!(;
        nc_dir=joinpath("data", "netcdf_stare"),
        lidarstemdir="./data",
        uv_path=joinpath("data/netcdf", "ekamsat_lidar_uv_20240428-20240613.nc"),
        nx=4000,
        nz=80,
        ntop=80,
        ic_list=nothing,
        overwrite=false,
        log_every=200)

    ctx = setup_sync_context_nc(; lidarstemdir=lidarstemdir, uv_path=uv_path, nc_dir=nc_dir, nx=nx, nz=nz)
    Env = ctx.Env
    Vn = ctx.Vn
    UV = ctx.UV
    beams = ctx.beams
    state = init_stream_state()

    ics = isnothing(ic_list) ? collect(eachindex(Env.ists)) : collect(ic_list)
    isempty(ics) && error("ic_list is empty")

    # Precreate mdv_sync variables once per file unless overwrite is requested.
    for ifile in eachindex(Env.ff)
        stem = splitext(basename(Env.ff[ifile]))[1]
        nc_path = joinpath(nc_dir, stem * ".nc")
        NCDataset(nc_path, "a") do ds
            if haskey(ds, "mdv_sync")
                if overwrite
                    ds["mdv_sync"][:] .= NaN
                end
            else
                v = defVar(ds, "mdv_sync", Float64, ("time",);
                    attrib=[
                        "long_name" => "sync-loop mdv from valid_chunk_mdv",
                        "units" => "m s-1",
                        "algorithm" => "LidarVNSync.extract_sync_window -> valid_chunk_mdv",
                        "noise_threshold" => NOISE_THR,
                        "ntop" => Int(ntop),
                        "_FillValue" => NaN,
                    ])
                v[:] .= NaN
            end
        end
    end

    n_total = length(ics)
    n_done = 0
    n_err = 0

    println("Pass-2 mdv_sync write: chunks=", n_total, ", nc_dir=", nc_dir)

    for (ii, ic) in enumerate(ics)
        try
            w = extract_sync_window(beams, Env, state, Vn, UV, ic; ntop=ntop, nc_dir=nc_dir)
            mdv_chunk = Float64.(w.mdv)

            # Write mdv across one or more files using global beam indices.
            g0 = w.ist
            while g0 <= w.ien
                ifile = findlast(Env.bigind_file_starts .<= g0)
                file_st = Env.bigind_file_starts[ifile]
                file_en = Env.bigind_file_ends[ifile]
                g1 = min(w.ien, file_en)

                local_st = g0 - file_st + 1
                local_en = g1 - file_st + 1
                chunk_st = g0 - w.ist + 1
                chunk_en = g1 - w.ist + 1

                stem = splitext(basename(Env.ff[ifile]))[1]
                nc_path = joinpath(nc_dir, stem * ".nc")
                NCDataset(nc_path, "a") do ds
                    ds["mdv_sync"][local_st:local_en] = mdv_chunk[chunk_st:chunk_en]
                end

                g0 = g1 + 1
            end

            n_done += 1
            if (ii % log_every == 0) || (ii == n_total)
                @printf("[%5d/%5d] done=%d err=%d\n", ii, n_total, n_done, n_err)
            end
        catch err
            n_err += 1
            @printf("ERROR chunk ic=%d: %s\n", ic, replace(sprint(showerror, err), '\n' => " | "))
        end
    end

    return (; nchunks=n_total, done=n_done, errors=n_err, ntop=ntop, nc_dir)
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

    sync = coarse_and_fine_lag(win.mdv, win.vn2_xcorr; prior_seconds=0.0)
    ref20 = refine_offset_20hz(win, sync, Vn)

    println("clock prior applied: ", win.lidar_clock_fast_by)
    @printf("coarse offset after iterative envelope passes = %.3f s\n", sync.coarse_offset)
    @printf("fine residual after coarse shift = %.3f s\n", sync.fine.lag_seconds)
    @printf("final offset (1 Hz stage) = %.3f s\n", sync.final_offset)
    @printf("20 Hz residual (raw) = %.3f s\n", ref20.final_offset_20hz)
    @printf("20 Hz residual (native step) = %.3f s [step=%.3f]\n", ref20.final_offset_native, ref20.native_step)
    @printf("final offset (1 Hz + 20 Hz residual) = %.3f s\n", sync.final_offset + ref20.final_offset_native)

    fig = figure(figsize=(10, 11))
    clf()

    subplot(5, 1, 1)
    pcolormesh((beams[:time][win.ist:win.ien] .% 1) .* 60, beams[:height] ./ 1e3, pd(win.beta), cmap=ColorMap("RdYlBu_r"))
    ylabel("height (km)")
    title("Chunk $(ic): backscatter")
    colorbar()

    subplot(5, 1, 2)
    plot(win.stare_dt, win.mdv, label="mdv")
    plot(win.stare_dt, win.vn2_xcorr, label="VelNED2 (1 Hz avg for sync)")
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
    ax1.plot(result.seq_ic, result.seq_final_native, marker="o", markersize=1,label="final offset (20 Hz native)")
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

"""
    write_daily_mdv_vn2!(; out_dir, nc_dir, lidarstemdir, uv_path, nx, nz, ntop,
                           ic_list, overwrite, log_path, log_every)

Production run: loops over all VN-covered lidar chunks, computes per-chunk sync
offsets (1 Hz sequential prior chain + 20 Hz native refinement), groups chunks by
the calendar date of their `stare_dt[1]`, and writes one NetCDF file per day to
`out_dir`.

Each daily file has **two** dimensions:

  time  (UNLIMITED) — one record per lidar beam; ncrcat-concatenable across daily files
  chunk (fixed)     — one record per sync chunk contributing to this calendar day

Variables under the `time` dimension:

  time            Int64   milliseconds since 1970-01-01T00:00:00 UTC
                          (CF-1.6; xarray/CDO auto-decode without keyword args;
                           Int64 ms avoids Float64 sub-ms rounding near the 2024 epoch)
  time_str        String  ISO-8601 UTC datetime of each beam
                          (e.g. "2024-04-29T23:59:01.020"; human-readable redundancy)
  mdv             Float32 SNR-weighted mean Doppler velocity                      [m s-1]
  vn2_aligned     Float32 VectorNav VelNED2 aligned to lidar time                 [m s-1]
  pitch           Float32 VectorNav pitch aligned to lidar time                   [degrees]
  roll            Float32 VectorNav roll aligned to lidar time                    [degrees]
  sync_chunk_id   Int32   which sync chunk this beam belongs to (one constant value
                          per chunk); cross-references chunk_id in the chunk dimension

Variables under the `chunk` dimension:

  chunk_id        Int32   global 1-based sync-chunk identifier (1 … ~5938 for the
                          full deployment); NOT a sample index within a chunk;
                          use to join to per-chunk diagnostics and to look up ist/ien
                          in the lidar index for the full Doppler velocity arrays
  lidar_t_start   Int64   epoch-ms of the first lidar beam in this sync chunk     [ms]
  lidar_t_end     Int64   epoch-ms of the last  lidar beam in this sync chunk     [ms]
  vn_t_start      Int64   epoch-ms of the first VectorNav 20 Hz sample in the chunk
                          window (fill=-1 when VN coverage was insufficient)      [ms]
  vn_t_end        Int64   epoch-ms of the last  VectorNav 20 Hz sample in the chunk
                          window (fill=-1 when VN coverage was insufficient)      [ms]
  offset_s        Float32 total VN timing offset applied to this chunk
                          (1 Hz sequential prior + 20 Hz native residual)         [s]

Progress and per-chunk quality are appended to `log_path` (CSV) after every chunk.
A summary line is printed to stdout every `log_every` chunks.
"""
function write_daily_mdv_vn2!(;
    out_dir      = joinpath("epsilon_data", "daily"),
    nc_dir       = "data/netcdf_stare",
    lidarstemdir = "./data",
    uv_path      = joinpath("data/netcdf", "ekamsat_lidar_uv_20240428-20240613.nc"),
    nx           = 4000,
    nz           = 80,
    ntop         = 80,
    ic_list      = nothing,
    overwrite    = false,
    log_path     = joinpath("epsilon_data",
                     "daily_mdv_vn2_$(Dates.format(Dates.now(), dateformat"yyyymmdd_HHMMSS")).log"),
    log_every    = 100,
)
    mkpath(out_dir)
    mkpath(dirname(log_path))

    ctx   = setup_sync_context_nc(; lidarstemdir, uv_path, nc_dir, nx, nz)
    Env   = ctx.Env
    Vn    = ctx.Vn
    UV    = ctx.UV
    beams = ctx.beams
    ics   = isnothing(ic_list) ? collect(ctx.icvn) : collect(ic_list)
    n_total = length(ics)
    @printf("write_daily_mdv_vn2!: %d chunks → %s\n  Log: %s\n", n_total, out_dir, log_path)

    # ── Phase 1: sequential 1 Hz offsets (establishes prior chain) ─────────
    println("Phase 1: 1 Hz sequential offsets …")
    t0 = time()
    seq_results = run_sequential_offsets(beams, Env, Vn, UV, ics; ntop=ntop, nc_dir=nc_dir)
    by_ic = Dict(r.ic => r for r in seq_results)
    @printf("  done in %.1f s\n", time() - t0)

    # ── Phase 2: 20 Hz refinement + beam-level data collection ─────────────
    println("Phase 2: 20 Hz refinement + data collection …")
    t0 = time()
    state = init_stream_state()

    # Pre-allocate per-chunk storage (NamedTuples with concrete arrays).
    chunk_data = NamedTuple[]
    sizehint!(chunk_data, n_total)

    open(log_path, "w") do logf
        println(logf, "ic,day,offset_1hz_s,offset_total_s,corr,vn_coverage,n_beams,nan_frac_vn2,status")
        flush(logf)

        for (ii, ic) in enumerate(ics)
            r      = get(by_ic, ic, nothing)
            status = "ok"
            corr   = NaN
            offset_total = NaN

            try
                w = extract_sync_window(beams, Env, state, Vn, UV, ic;
                        ntop=ntop, nc_dir=nc_dir)
                day     = Date(w.stare_dt[1])
                n_beams = length(w.stare_dt)
                dt1     = dt_seconds(w.stare_dt)
                (!isfinite(dt1) || dt1 <= 0) && (dt1 = TIMESTEP)

                offset_1hz = isnothing(r) ? NaN :
                             (isfinite(r.final_offset) ? r.final_offset : NaN)
                bad_sync = (w.vn_coverage < 0.5) || isnothing(r) ||
                           !get(r, :accepted_sync, true) || !isfinite(offset_1hz)

                if bad_sync
                    status = w.vn_coverage < 0.5 ? "low_vn_cov" :
                             isnothing(r)         ? "no_result"  : "rejected_sync"
                    push!(chunk_data, (;
                        ic, day,
                        lidar_dt    = w.stare_dt,
                        mdv         = Float32.(w.mdv),
                        vn2_aligned = fill(Float32(NaN), n_beams),
                        pitch       = Float32.(w.pitch),
                        roll        = Float32.(w.roll),
                        offset_s    = Float32(NaN),
                        vn_t_start  = isempty(w.vndt20_chunk) ? missing : w.vndt20_chunk[1],
                        vn_t_end    = isempty(w.vndt20_chunk) ? missing : w.vndt20_chunk[end],
                    ))
                else
                    # Reuse the 1 Hz result from Phase 1; only run 20 Hz here.
                    s1_stub = (; prior_seconds = isfinite(r.prior_offset) ? r.prior_offset : 0.0,
                                 final_offset  = offset_1hz)
                    rf = refine_offset_20hz(w, s1_stub, Vn)
                    resid        = isfinite(rf.final_offset_native) ? rf.final_offset_native : 0.0
                    offset_total = offset_1hz + resid
                    vn2_aligned  = Float32.(rf.vn2_1s_aligned)
                    corr         = finite_overlap_corr(w.mdv, vn2_aligned)
                    status       = isfinite(corr) && corr > 0.3 ? "ok" : "low_corr"

                    push!(chunk_data, (;
                        ic, day,
                        lidar_dt    = w.stare_dt,
                        mdv         = Float32.(w.mdv),
                        vn2_aligned,
                        pitch       = Float32.(rf.pitch_1s_aligned),
                        roll        = Float32.(rf.roll_1s_aligned),
                        offset_s    = Float32(offset_total),
                        vn_t_start  = isempty(w.vndt20_chunk) ? missing : w.vndt20_chunk[1],
                        vn_t_end    = isempty(w.vndt20_chunk) ? missing : w.vndt20_chunk[end],
                    ))
                end

                @printf(logf, "%d,%s,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%s\n",
                    ic, string(day),
                    isfinite(offset_1hz)   ? offset_1hz   : -9999.0,
                    isfinite(offset_total) ? offset_total : -9999.0,
                    isfinite(corr)         ? corr         : -9999.0,
                    w.vn_coverage, n_beams, w.vn2_xcorr_nan_frac, status)
                flush(logf)

                if ii % log_every == 0 || ii == n_total
                    @printf("  [%5d/%5d] ic=%-5d %-12s off=%7.3f corr=%5.3f  %s\n",
                        ii, n_total, ic, string(day),
                        isfinite(offset_total) ? offset_total : -9999.0,
                        isfinite(corr) ? corr : -9999.0, status)
                end
            catch err
                status = "error:" * replace(sprint(showerror, err), '\n' => " | ")
                @printf(logf, "%d,---,-9999,-9999,-9999,-9999,0,-9999,%s\n", ic, status)
                flush(logf)
                @printf("ERROR ic=%d: %s\n", ic, sprint(showerror, err))
            end
        end
    end
    @printf("  done in %.1f s\n", time() - t0)

    # ── Phase 3: write one NC file per calendar day ─────────────────────────
    println("Phase 3: writing daily NetCDF files …")
    days   = sort(unique([c.day for c in chunk_data]))
    n_days = length(days)
    n_ok   = count(c -> any(isfinite, c.vn2_aligned), chunk_data)
    @printf("  %d chunks → %d days  (%d chunks with valid vn2_aligned)\n",
        length(chunk_data), n_days, n_ok)

    for day in days
        nc_out = joinpath(out_dir,
            "mdv_vn2_sync_$(Dates.format(day, dateformat"yyyymmdd")).nc")
        if isfile(nc_out) && !overwrite
            @printf("  Skip (exists): %s\n", basename(nc_out))
            continue
        end

        day_chunks   = [c for c in chunk_data if c.day == day]
        n_chunks_day = length(day_chunks)

        # ── per-beam (time dimension) arrays ─────────────────────────────────
        lidar_dt_all = vcat([c.lidar_dt    for c in day_chunks]...)
        mdv_all      = vcat([c.mdv         for c in day_chunks]...)
        vn2_all      = vcat([c.vn2_aligned for c in day_chunks]...)
        pitch_all    = vcat([c.pitch        for c in day_chunks]...)
        roll_all     = vcat([c.roll         for c in day_chunks]...)
        scid_all     = vcat([fill(Int32(c.ic), length(c.lidar_dt)) for c in day_chunks]...)
        nrec         = length(lidar_dt_all)

        epoch_ref_ms  = Dates.datetime2epochms(DateTime(1970, 1, 1, 0, 0, 0))
        FILL_MS       = Int64(-1)   # impossible for valid Unix ms timestamps (all > 0)
        to_ms(t::DateTime) = Int64(Dates.datetime2epochms(t) - epoch_ref_ms)
        to_ms(::Missing)   = FILL_MS

        time_ms  = [to_ms(t) for t in lidar_dt_all]
        time_str = [Dates.format(t, dateformat"yyyy-mm-ddTHH:MM:SS.sss") for t in lidar_dt_all]

        # ── per-chunk (chunk dimension) arrays ───────────────────────────────
        chunk_ids    = Int32[c.ic                        for c in day_chunks]
        lidar_t_s_ms = [to_ms(c.lidar_dt[1])            for c in day_chunks]
        lidar_t_e_ms = [to_ms(c.lidar_dt[end])          for c in day_chunks]
        vn_t_s_ms    = [to_ms(c.vn_t_start)             for c in day_chunks]
        vn_t_e_ms    = [to_ms(c.vn_t_end)               for c in day_chunks]
        off_chunk    = Float32[c.offset_s                for c in day_chunks]

        # ── write NetCDF ──────────────────────────────────────────────────────
        isfile(nc_out) && rm(nc_out)
        ds = NCDataset(nc_out, "c")
        defDim(ds, "time",  Inf)            # UNLIMITED → ncrcat-concatenable across days
        defDim(ds, "chunk", n_chunks_day)   # one record per sync chunk in this day

        ms_units = "milliseconds since 1970-01-01T00:00:00 UTC"
        fv32     = Float32(NaN)

        # time-dimension variables
        v_t      = defVar(ds, "time",          Int64,   ("time",);
            attrib=["units"=>ms_units, "calendar"=>"standard",
                    "standard_name"=>"time", "axis"=>"T",
                    "long_name"=>"milliseconds since 1970-01-01T00:00:00 UTC"])
        v_tstr   = defVar(ds, "time_str",      String,  ("time",);
            attrib=["long_name"=>"ISO-8601 UTC datetime of each lidar beam (millisecond precision)"])
        v_mdv    = defVar(ds, "mdv",           Float32, ("time",);
            attrib=["units"=>"m s-1",   "long_name"=>"SNR-weighted mean Doppler velocity",
                    "_FillValue"=>fv32])
        v_vn2    = defVar(ds, "vn2_aligned",   Float32, ("time",);
            attrib=["units"=>"m s-1",   "long_name"=>"VectorNav VelNED2 aligned to lidar time",
                    "_FillValue"=>fv32])
        v_pit    = defVar(ds, "pitch",         Float32, ("time",);
            attrib=["units"=>"degrees", "long_name"=>"VectorNav pitch aligned to lidar time",
                    "_FillValue"=>fv32])
        v_rol    = defVar(ds, "roll",          Float32, ("time",);
            attrib=["units"=>"degrees", "long_name"=>"VectorNav roll aligned to lidar time",
                    "_FillValue"=>fv32])
        v_scid   = defVar(ds, "sync_chunk_id", Int32,   ("time",);
            attrib=["long_name"=>"which sync chunk this beam belongs to (one constant value per chunk); cross-references chunk_id in the chunk dimension"])

        # chunk-dimension variables
        v_cid  = defVar(ds, "chunk_id",      Int32,  ("chunk",);
            attrib=["long_name"=>"global 1-based sync-chunk identifier (1…~5938 for the full deployment); NOT a sample index within a chunk; join on sync_chunk_id to link beams to chunks"])
        v_lst  = defVar(ds, "lidar_t_start", Int64,  ("chunk",);
            attrib=["units"=>ms_units, "calendar"=>"standard",
                    "long_name"=>"epoch-ms of the first lidar beam in this sync chunk"])
        v_len  = defVar(ds, "lidar_t_end",   Int64,  ("chunk",);
            attrib=["units"=>ms_units, "calendar"=>"standard",
                    "long_name"=>"epoch-ms of the last lidar beam in this sync chunk"])
        v_vst  = defVar(ds, "vn_t_start",    Int64,  ("chunk",);
            attrib=["units"=>ms_units, "calendar"=>"standard", "_FillValue"=>FILL_MS,
                    "long_name"=>"epoch-ms of the first VectorNav 20 Hz sample in the chunk window (fill=-1 when VN coverage insufficient)"])
        v_ven  = defVar(ds, "vn_t_end",      Int64,  ("chunk",);
            attrib=["units"=>ms_units, "calendar"=>"standard", "_FillValue"=>FILL_MS,
                    "long_name"=>"epoch-ms of the last VectorNav 20 Hz sample in the chunk window (fill=-1 when VN coverage insufficient)"])
        v_off  = defVar(ds, "offset_s",      Float32, ("chunk",);
            attrib=["units"=>"s", "_FillValue"=>fv32,
                    "long_name"=>"total VN timing offset applied to this chunk (1 Hz sequential prior + 20 Hz native residual)"])

        # write time dimension
        v_t[:]     = time_ms
        v_tstr[:]  = time_str
        v_mdv[:]   = mdv_all
        v_vn2[:]   = vn2_all
        v_pit[:]   = pitch_all
        v_rol[:]   = roll_all
        v_scid[:]  = scid_all

        # write chunk dimension
        v_cid[:] = chunk_ids
        v_lst[:] = lidar_t_s_ms
        v_len[:] = lidar_t_e_ms
        v_vst[:] = vn_t_s_ms
        v_ven[:] = vn_t_e_ms
        v_off[:] = off_chunk

        ds.attrib["date"]        = string(day)
        ds.attrib["n_chunks"]    = n_chunks_day
        ds.attrib["n_beams"]     = nrec
        ds.attrib["ntop"]        = Int(ntop)
        ds.attrib["nc_dir"]      = nc_dir
        ds.attrib["created_utc"] = string(Dates.now(Dates.UTC))
        ds.attrib["Conventions"] = "CF-1.6"
        close(ds)

        @printf("  Wrote: %-52s  chunks=%d  beams=%d\n",
            basename(nc_out), length(day_chunks), nrec)
    end

    return (; n_total, n_days, out_dir, log_path)
end

function diagnostic_example_chunks(beams, Env, Vn, UV, ic_list; ntop=80, nc_dir=nothing)
    state = init_stream_state()
    ics = collect(ic_list)
    isempty(ics) && error("ic_list is empty")

    nplot = length(ics)
    fig = figure(figsize=(11, 2.6 * nplot))
    clf()
    rows = max(nplot, 1)
    results = NamedTuple[]

    for (i, ic) in enumerate(ics)
        win = extract_sync_window(beams, Env, state, Vn, UV, ic; ntop=ntop, nc_dir=nc_dir)
        sync = coarse_and_fine_lag(win.mdv, win.vn2_xcorr; prior_seconds=0.0)
        ref20 = refine_offset_20hz(win, sync, Vn)
        corr_aligned = finite_overlap_corr(win.mdv, ref20.vn2_1s_aligned)

        subplot(rows, 1, i)
        plot(win.stare_dt, win.mdv, label="mdv", linewidth=1.6)
        plot(win.stare_dt, win.vn2_xcorr, label="vn2 raw 1 s", linewidth=1.0, alpha=0.6)
        plot(win.stare_dt, ref20.vn2_1s_aligned, label="vn2 aligned 1 s", linewidth=1.4)
        plot(win.stare_dt, ref20.mdv_residual_1s, label="mdv - vn2 aligned", linewidth=1.0, alpha=0.75)
        ylabel("m s^-1")
        title(@sprintf("chunk %d | 1 Hz offset = %.3f s | 20 Hz residual = %.3f s | total = %.3f s | corr = %.3f", ic, sync.final_offset, ref20.final_offset_20hz, sync.final_offset + ref20.final_offset_native, corr_aligned))
        grid(true)
        i == 1 && legend(loc="best")
        i == rows && xlabel("time")

        push!(results, (; ic, win, sync, ref20, corr_aligned))
    end

    tight_layout()
    return (; fig, results)
end

end