using Pkg
Pkg.activate(".")

using Dates
using JLD2
using Statistics

const EPS_DIR = "epsilon_data"
const OUT_DIR = joinpath(EPS_DIR, "hmix")

"Height at lidar range gate center in meters above the lidar transmitter."
rng_height(iz::Integer; rangegate::Float64 = 24.0) = rangegate * (iz - 1 + 0.5)

"Return true only for physically valid epsilon values."
is_valid_eps(x) = !ismissing(x) && isfinite(x) && (x > 0)

"Compute first-crossing h_mix index and height for one epsilon profile."
function hmix_first_crossing(profile::AbstractVector, height::AbstractVector; min_valid_below::Int = 1)
    nz = length(profile)
    nz == length(height) || throw(ArgumentError("profile and height lengths must match"))

    h_idx = missing
    h_m = missing

    for iz in 2:nz
        eps_i = profile[iz]
        if !is_valid_eps(eps_i)
            continue
        end

        below_vals = [profile[j] for j in 1:(iz - 1) if is_valid_eps(profile[j])]
        if length(below_vals) < min_valid_below
            continue
        end

        thr = mean(below_vals) / 3
        if eps_i <= thr
            h_idx = iz
            h_m = height[iz]
            break
        end
    end

    return h_idx, h_m
end

"Extract epsilon matrix and metadata from either supported file schema."
function extract_epsilon_and_meta(path::AbstractString)
    d = load(path)

    eps_key = if haskey(d, "epsilon")
        "epsilon"
    elseif haskey(d, "epsi")
        "epsi"
    else
        error("No epsilon matrix key found in $path")
    end

    epsilon = d[eps_key]
    ndims(epsilon) == 2 || error("Expected 2D epsilon matrix in $path")

    nprof, nz = size(epsilon)
    height = if haskey(d, "height")
        h = d["height"]
        length(h) == nz || error("height length mismatch in $path")
        collect(h)
    else
        [rng_height(iz) for iz in 1:nz]
    end

    times_start = if haskey(d, "start_dt")
        d["start_dt"]
    elseif haskey(d, "dtime_st")
        d["dtime_st"]
    else
        nothing
    end

    times_end = if haskey(d, "end_dt")
        d["end_dt"]
    elseif haskey(d, "dtime_en")
        d["dtime_en"]
    else
        nothing
    end

    return epsilon, height, times_start, times_end, nprof
end

"Compute h_mix arrays for every profile row in epsilon matrix."
function compute_hmix_for_matrix(epsilon::AbstractMatrix, height::AbstractVector; min_valid_below::Int = 1)
    nprof = size(epsilon, 1)

    hmix_index = Vector{Union{Missing, Int}}(missing, nprof)
    hmix_height_m = Vector{Union{Missing, Float64}}(missing, nprof)
    n_valid_eps = zeros(Int, nprof)

    for ip in 1:nprof
        prof = vec(epsilon[ip, :])
        n_valid_eps[ip] = count(is_valid_eps, prof)
        idx, hm = hmix_first_crossing(prof, height; min_valid_below = min_valid_below)
        hmix_index[ip] = idx
        hmix_height_m[ip] = hm
    end

    return hmix_index, hmix_height_m, n_valid_eps
end

"Process one epsilon file and write h_mix output file."
function process_file(path::AbstractString; out_dir::AbstractString = OUT_DIR, min_valid_below::Int = 1)
    epsilon, height, tstart, tend, nprof = extract_epsilon_and_meta(path)
    hmix_index, hmix_height_m, n_valid_eps = compute_hmix_for_matrix(epsilon, height; min_valid_below = min_valid_below)

    mkpath(out_dir)
    base = splitext(basename(path))[1]
    out_path = joinpath(out_dir, "hmix_$(base).jld2")

    n_found = count(!ismissing, hmix_index)
    n_missing = nprof - n_found
    med_height = begin
        vals = collect(skipmissing(hmix_height_m))
        isempty(vals) ? missing : median(vals)
    end

    out = Dict{String, Any}(
        "source_file" => path,
        "computed_at" => now(),
        "criterion" => "first iz where epsilon[iz] <= mean(epsilon[1:iz-1] valid_only)/3",
        "min_valid_below" => min_valid_below,
        "profile_count" => nprof,
        "height_m" => height,
        "hmix_index" => hmix_index,
        "hmix_height_m" => hmix_height_m,
        "n_valid_eps" => n_valid_eps,
        "n_hmix_found" => n_found,
        "n_hmix_missing" => n_missing,
        "median_hmix_height_m" => med_height,
    )

    if tstart !== nothing
        out["time_start"] = tstart
    end
    if tend !== nothing
        out["time_end"] = tend
    end

    save(out_path, out)

    println("processed: $(basename(path))")
    println("  profiles: $nprof, h_mix found: $n_found, missing: $n_missing")
    println("  output: $out_path")
    return out_path
end

"Return source epsilon files to process from epsilon_data directory."
function list_epsilon_files(eps_dir::AbstractString = EPS_DIR)
    files = readdir(eps_dir)
    selected = filter(files) do f
        (startswith(f, "epsilon_") && endswith(f, ".jld2")) ||
        (startswith(f, "epsi_stare_chunks_") && endswith(f, ".jld"))
    end
    sort!(selected)
    return [joinpath(eps_dir, f) for f in selected]
end

function main(; min_valid_below::Int = 1)
    files = list_epsilon_files()
    println("found $(length(files)) epsilon files")

    out_paths = String[]
    for path in files
        push!(out_paths, process_file(path; min_valid_below = min_valid_below))
    end

    println("done. wrote $(length(out_paths)) h_mix files to $(OUT_DIR)")
end

main()