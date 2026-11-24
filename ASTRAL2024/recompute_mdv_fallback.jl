using Pkg
Pkg.activate(@__DIR__)

using Dates
using Statistics
using NCDatasets
using JLD2

include(joinpath(@__DIR__, "DopplerTurbulence.jl"))
include(joinpath(@__DIR__, "lidar_vn_sync.jl"))
using .DopplerTurbulence

good(x) = !ismissing(x) && isfinite(x)
tonan(x) = good(x) ? Float64(x) : NaN

function interpolate_at(times, values, t)
    j = searchsortedfirst(times, t)
    if j <= length(times) && times[j] == t
        return values[j]
    elseif j <= 1 || j > length(times)
        return NaN
    end
    a, b = values[j - 1], values[j]
    good(a) && good(b) || return NaN
    f = Dates.value(t - times[j - 1]) / Dates.value(times[j] - times[j - 1])
    return (1 - f) * a + f * b
end

function load_wind_context(sonde_path, ship_path)
    sonde = NCDataset(sonde_path) do ds
        times = ds["time"][:]
        u = tonan.(ds["u"][:, :])
        v = tonan.(ds["v"][:, :])
        heights = tonan.(ds["gp_height"][:])
        bad = findall(vec(any(>(25), u; dims=1) .| any(>(25), v; dims=1)))
        length(bad) == 2 || error("Expected the two production sounding QC cases; found $(length(bad))")
        for (n, i) in enumerate(bad)
            1 < i < length(times) || error("Cannot repair sounding at index $i")
            levels = n == 1 ? (10:size(u, 1)) : (1:size(u, 1))
            u[levels, i] .= (u[levels, i - 1] .+ u[levels, i + 1]) ./ 2
            v[levels, i] .= (v[levels, i - 1] .+ v[levels, i + 1]) ./ 2
        end
        (; times, u, v, heights)
    end
    ship = NCDataset(ship_path) do ds
        times = ds["time"][:]
        sog = tonan.(ds["sog"][:])
        cog = tonan.(ds["cog"][:])
        u = [good(s) && good(c) && s <= 20 ? s * sind(c) : NaN for (s, c) in zip(sog, cog)]
        v = [good(s) && good(c) && s <= 20 ? s * cosd(c) : NaN for (s, c) in zip(sog, cog)]
        (; times, u, v)
    end
    issorted(sonde.times) && issorted(ship.times) || error("Wind timestamps are not sorted")
    return (; sonde, ship, daily=Dict{Date, Any}())
end

function daily_wind!(context, day)
    get!(context.daily, day) do
        path = joinpath(@__DIR__, "data", "netcdf", "v1",
                        "ekamsat_lidar_uv_$(Dates.format(day, dateformat"yyyymmdd"))_v1.nc")
        if !isfile(path)
            @warn "No daily lidar-wind file; using sounding winds where available" day path
            return nothing
        end
        NCDataset(path) do ds
            times = ds["time"][:]
            ur, vr = tonan.(ds["ur"][:, :]), tonan.(ds["vr"][:, :])
            issorted(times) && !isempty(times) && all(t -> Date(t) == day, times) ||
                error("Unexpected time coverage in $path")
            size(ur) == size(vr) && size(ur, 1) == length(times) ||
                error("Unexpected wind dimensions in $path")
            (; times, ur, vr)
        end
    end
end

blend(a, b) = good(a) && good(b) ? (a + b) / 2 : good(a) ? a : good(b) ? b : NaN

function chunk_winds(context, times, nz)
    t = first(times)
    ship_u = interpolate_at(context.ship.times, context.ship.u, t)
    ship_v = interpolate_at(context.ship.times, context.ship.v, t)
    sonde_u = [interpolate_at(context.sonde.times, view(context.sonde.u, iz, :), t) - ship_u for iz in 1:nz]
    sonde_v = [interpolate_at(context.sonde.times, view(context.sonde.v, iz, :), t) - ship_v for iz in 1:nz]
    ur = fill(NaN, length(times), nz)
    vr = fill(NaN, length(times), nz)
    for (i, time) in enumerate(times)
        daily = daily_wind!(context, Date(time))
        if !isnothing(daily)
            size(daily.ur, 2) >= nz || error("Daily lidar wind has fewer than $nz levels")
            j = min(searchsortedfirst(daily.times, time), length(daily.times))
            for iz in 1:nz
                ur[i, iz] = daily.ur[j, iz]
                vr[i, iz] = daily.vr[j, iz]
            end
        end
        for iz in 1:nz
            ur[i, iz] = blend(ur[i, iz], sonde_u[iz])
            vr[i, iz] = blend(vr[i, iz], sonde_v[iz])
        end
    end
    return ur, vr
end

function aligned_vn(ic, times, n)
    path = joinpath(@__DIR__, "data", "vn_sync_chunk_daily",
                    "mdv_vn2_sync_$(Dates.format(Date(first(times)), dateformat"yyyymmdd")).nc")
    isfile(path) || return nothing
    return NCDataset(path) do ds
        k = findfirst(==(ic), ds["chunk_ic"][:])
        isnothing(k) && return nothing
        a, b = Int(ds["time_start_idx"][k]), Int(ds["time_end_idx"][k])
        b - a + 1 == n || error("VectorNav beam count differs for chunk $ic")
        (; heave=tonan.(ds["vn2_aligned"][a:b]),
           pitch=tonan.(ds["pitch"][a:b]) .* (pi / 180),
           roll=tonan.(ds["roll"][a:b]) .* (pi / 180))
    end
end

function fit_epsilon(w, pitch, roll, ur, vr)
    sigma = 2pi / 180
    result = D2_rho_small_disp(w, pitch, roll, ur, vr, sigma, sigma;
                               Lmin=48.0, Lmax=125.0, min_fit_bins=3)
    slopes = result[3]
    epsi = [ismissing(a) || !isfinite(a) ? -3.0 : a < 0 ? -1.0 : epsilon(a) for a in slopes]
    for iz in eachindex(epsi)
        any(i -> good(ur[i, iz]) && good(vr[i, iz]), axes(ur, 1)) || (epsi[iz] = -4.0)
    end
    return (; epsi, A=slopes, nbins=result[7], nrho=result[13],
            rho=result[2], D2=result[1])
end

function paired_variance(a, b)
    va, vb = fill(NaN, size(a, 2)), fill(NaN, size(a, 2))
    count_paired = zeros(Int, size(a, 2))
    for iz in axes(a, 2)
        ii = [i for i in axes(a, 1) if good(a[i, iz]) && good(b[i, iz])]
        count_paired[iz] = length(ii)
        if length(ii) >= 2
            va[iz], vb[iz] = var(a[ii, iz]), var(b[ii, iz])
        end
    end
    return va, vb, count_paired
end

function main(args)
    4 <= length(args) <= 6 || error(
        "Usage: julia --project=. recompute_mdv_fallback.jl START_DATE END_DATE_EXCLUSIVE SHIP_NAV_NC OUTPUT_JLD2 [CHUNK_ID] [--compare-vn]")
    first_day, end_day = Date(args[1]), Date(args[2])
    first_day < end_day || error("END_DATE_EXCLUSIVE must follow START_DATE")
    ship_path, outfile = abspath(args[3]), abspath(args[4])
    isfile(ship_path) || error("Ship navigation file not found: $ship_path")
    ispath(outfile) && error("Output already exists: $outfile")
    rest = args[5:end]
    compare_vn = "--compare-vn" in rest
    chunk_args = filter(!=("--compare-vn"), rest)
    length(chunk_args) <= 1 || error("At most one chunk ID may be specified")

    env = LidarVNSync.load_lidar_indices_and_files(; lidarstemdir=joinpath(@__DIR__, "data"))
    chunk_ids = findall(i -> first_day <= Date(env.dtime[env.ists[i]]) < end_day, eachindex(env.ists))
    if !isempty(chunk_args)
        ic = parse(Int, only(chunk_args))
        ic in chunk_ids || error("Chunk $ic is outside the requested date range")
        chunk_ids = [ic]
    end
    isempty(chunk_ids) && error("No lidar chunks in the requested date range")
    context = load_wind_context(joinpath(@__DIR__, "data", "netcdf", "EKAMSAT_RAOB_TGT_combined_80lev.nc"), ship_path)
    nz, n = 80, length(chunk_ids)
    beams = LidarVNSync.init_periodic_beams(4000, nz)
    state = LidarVNSync.init_stream_state()
    epsi = fill(-3.0, n, nz)
    slopes = Matrix{Union{Missing,Float64}}(missing, n, nz)
    nbins, nrho = zeros(Int, n, nz), zeros(Int, n, nz)
    rho = Array{Union{Missing,Float64}}(missing, n, 100, nz)
    D2 = similar(rho)
    fill!(D2, missing)
    gps_start = Vector{DateTime}(undef, n)
    variance_mdv_all = fill(NaN, n, nz)
    mdv_samples = zeros(Int, n, nz)
    variance_mdv = fill(NaN, n, nz)
    variance_vn = fill(NaN, n, nz)
    paired_samples = zeros(Int, n, nz)
    epsi_mdv_paired = fill(NaN, n, nz)
    epsi_vn_paired = fill(NaN, n, nz)
    height = Float64[]

    for (row, ic) in enumerate(chunk_ids)
        info = LidarVNSync.ensure_chunk_loaded_nc!(beams, env, state, ic)
        raw_times = env.dtime[info.ist:info.ien]
        clock_offset = Millisecond(round(Int, 1000 * LidarVNSync.chunks.fit_offset(first(raw_times))))
        times = raw_times .- clock_offset
        gps_start[row] = first(times)
        height = tonan.(beams[:height][:])
        all(isapprox.(height, context.sonde.heights; atol=0.01)) ||
            error("Sounding and lidar height grids differ")
        dv_raw = beams[:dopplervel][info.ist:info.ien, 1:nz]
        intensity = beams[:intensity][info.ist:info.ien, 1:nz]
        dv = [good(v) && good(s) && s >= LidarVNSync.NOISE_THR ? Float64(v) : NaN
              for (v, s) in zip(dv_raw, intensity)]
        dv = reshape(dv, size(dv_raw))
        levels = findall(vec(all(isfinite.(dv); dims=1)))
        if isempty(levels)
            @warn "No complete gates for MDV" ic
            continue
        end
        heave = vec(mean(dv[:, levels]; dims=2))
        w = dv .- heave
        for iz in 1:nz
            valid = filter(isfinite, w[:, iz])
            mdv_samples[row, iz] = length(valid)
            if length(valid) >= 2
                variance_mdv_all[row, iz] = var(valid)
            end
        end
        ur, vr = chunk_winds(context, times, nz)
        fit = fit_epsilon(w, zeros(length(times)), zeros(length(times)), ur, vr)
        epsi[row, :] .= fit.epsi
        slopes[row, :] .= fit.A
        nbins[row, :] .= fit.nbins
        nrho[row, :] .= fit.nrho
        rho[row, :, :] .= fit.rho
        D2[row, :, :] .= fit.D2
        if compare_vn
            vn = aligned_vn(ic, times, length(times))
            if isnothing(vn)
                @warn "No aligned VectorNav data for comparison" ic
            else
                ok = isfinite.(vn.heave) .& isfinite.(vn.pitch) .& isfinite.(vn.roll)
                w_mdv = copy(w)
                w_mdv[.!ok, :] .= NaN
                w_vn = dv .- vn.heave
                w_vn[.!ok, :] .= NaN
                variance_mdv[row, :], variance_vn[row, :], paired_samples[row, :] = paired_variance(w_mdv, w_vn)
                if count(ok) >= 60
                    p = ifelse.(ok, vn.pitch, 0.0)
                    r = ifelse.(ok, vn.roll, 0.0)
                    epsi_mdv_paired[row, :] .= fit_epsilon(w_mdv, zeros(length(times)), zeros(length(times)), ur, vr).epsi
                    epsi_vn_paired[row, :] .= fit_epsilon(w_vn, p, r, ur, vr).epsi
                else
                    @warn "Too few aligned VectorNav beams for dissipation comparison" ic available=count(ok)
                end
            end
        end
        println("chunk $ic $(first(raw_times)): $(count(>(0), fit.epsi)) positive gates")
        flush(stdout)
    end
    save(outfile, Dict(
        "epsi" => epsi, "A" => slopes, "nbins" => nbins, "nrho" => nrho,
        "rho" => rho, "D2" => D2, "height" => height, "chunk_id" => chunk_ids,
        "dtime_st" => env.dtime[env.ists[chunk_ids]],
        "dtime_en" => env.dtime[env.iens[chunk_ids]], "dtime_st_gps" => gps_start,
        "variance_mdv" => variance_mdv_all, "mdv_samples" => mdv_samples,
        "variance_mdv_paired" => variance_mdv, "variance_vn_paired" => variance_vn,
        "paired_samples" => paired_samples,
        "epsi_mdv_paired" => epsi_mdv_paired, "epsi_vn_paired" => epsi_vn_paired,
        "wind_method" => "equal lidar and ship-relative sonde when both valid; otherwise available source",
        "lidar_wind_height_method" => "production gate-index matching",
        "mdv_method" => "mean Doppler velocity over gates complete for the chunk",
        "min_separation_m" => 48.0, "max_separation_m" => 125.0,
        "pitch_roll" => "zero for MDV; aligned VectorNav for comparison",
        "ship_navigation_file" => ship_path,
    ))
    println("Wrote $outfile")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
