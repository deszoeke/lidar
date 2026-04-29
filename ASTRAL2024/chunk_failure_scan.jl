"""
chunk_failure_scan.jl

Chunk-by-chunk failure diagnosis for lidar–VectorNav sync, starting at chunk 280.
For each chunk in the scan range, prints and saves:

  ifile          - lidar file index for this chunk
  new_file       - true when this chunk is in a different file than the previous one
  vn_coverage    - fraction of the chunk's time span covered by VN 20 Hz data
                   (< 0.5  → sync is skipped upstream by run_sequential_offsets)
  vn_nan_frac    - fraction of vn2_xcorr values that are NaN
                   (> 0.15 → sync is skipped upstream)
  ngoodlevels    - number of range gates with clean MDV (0 → no MDV signal)
  mdv_nfinite    - number of finite MDV time samples available for xcorr
  lidar_dt_s     - median beam-to-beam spacing in this chunk [s]
  lidar_dt_std_s - std-dev of beam spacing (large → lidar clock discontinuity)
  vn_dt_s        - median VN sample spacing inside the chunk window [s]
  vn_dt_std_s    - std-dev of VN spacing (large → VN clock jitter / gap)
  vn_gap_max_s   - largest single VN gap [s]  (> 2 s is significant)
  offset_s       - timing offset estimate [s]  (-9999 = sentinel / failed)
  corr           - Pearson corr(mdv, vn2_shifted) after applying offset
  coh_5_20s      - Welch coherence, 5-20 s band, after applying offset
  flag           - pipe-separated failure tags (ok if nothing flagged)

Usage from the project root:
  julia --project=. chunk_failure_scan.jl
"""

using Pkg
Pkg.activate(".")

using Dates
using Printf
using Statistics
using FFTW
using NCDatasets

include("./lidar_vn_sync.jl")
using .LidarVNSync

# ── helpers (mirror of workbench notebook) ─────────────────────────────────

function pearson_finite(x, y)
    xv, yv, _ = LidarVNSync.finite_overlap(Float64.(x), Float64.(y))
    n = length(xv)
    n < 4 && return NaN
    sx = std(xv); sy = std(yv)
    (!isfinite(sx) || !isfinite(sy) || sx == 0 || sy == 0) && return NaN
    return cor(xv, yv)
end

function band_coherence_welch_5_20s(x, y; dt=1.0, flo=1/20, fhi=1/5)
    xv, yv, _ = LidarVNSync.finite_overlap(Float64.(x), Float64.(y))
    n = length(xv)
    n < 48 && return NaN
    nseg = min(128, n); nseg < 32 && return NaN
    step = max(8, nseg ÷ 2)
    win  = 0.5 .- 0.5 .* cos.(2pi .* (0:(nseg-1)) ./ (nseg-1))
    nfreq = div(nseg, 2) + 1
    Sxx = zeros(Float64, nfreq); Syy = zeros(Float64, nfreq)
    Sxy = zeros(ComplexF64, nfreq); k = 0
    for i0 in 1:step:(n-nseg+1)
        xs = xv[i0:(i0+nseg-1)] .- mean(xv[i0:(i0+nseg-1)])
        ys = yv[i0:(i0+nseg-1)] .- mean(yv[i0:(i0+nseg-1)])
        X = FFTW.rfft(xs .* win); Y = FFTW.rfft(ys .* win)
        Sxx .+= abs2.(X); Syy .+= abs2.(Y); Sxy .+= X .* conj.(Y); k += 1
    end
    k == 0 && return NaN
    Sxx ./= k; Syy ./= k; Sxy ./= k
    coh  = abs2.(Sxy) ./ (Sxx .* Syy .+ eps())
    fs   = 1/dt
    freq = collect(0:div(nseg,2)) .* (fs/nseg)
    band = (freq .>= flo) .& (freq .<= fhi)
    any(band) || return NaN
    c = filter(isfinite, coh[band])
    isempty(c) && return NaN
    return mean(c)
end

# ── setup ──────────────────────────────────────────────────────────────────

ctx   = LidarVNSync.setup_sync_context()
Env   = ctx.Env
Vn    = ctx.Vn
UV    = ctx.UV
beams = ctx.beams
nz    = ctx.nz

# ── scan parameters ────────────────────────────────────────────────────────

SCAN_START    = 278          # start a couple before 280 for context
SCAN_END      = 310          # extend past the first sentinel (chunk 300)
SYNC_CORR_THR = 0.5
SYNC_COH_THR  = 0.9
NC_DIR        = "./data/netcdf_stare"

scan_chunks = SCAN_START:SCAN_END

# ── run 1-Hz sync over the full scan range so the prior chain is intact ────

println("Running 1-Hz sync for chunks $(SCAN_START):$(SCAN_END) …")
res = LidarVNSync.process_sync_data(
    beams, Env, Vn, UV, collect(scan_chunks);
    ntop=nz,
    reset_nan_log=false,
)
offset_by_ic = Dict(zip(res.seq_ic, Float64.(res.seq_final_native)))

# ── per-chunk diagnostic loop ──────────────────────────────────────────────

header = join(["chunk","ifile","new_file","vn_coverage","vn_nan_frac",
               "ngoodlevels","mdv_nfinite",
               "lidar_dt_s","lidar_dt_std_s",
               "vn_dt_s","vn_dt_std_s","vn_gap_max_s",
               "offset_s","corr","coh_5_20s","flag"], ",")
println("\n" * header)

state     = LidarVNSync.init_stream_state()
prev_ifile = -1
rows = NamedTuple[]

for ic in scan_chunks

    win = LidarVNSync.extract_sync_window(
        beams, Env, state, Vn, UV, ic;
        ntop=nz, nc_dir=NC_DIR)

    # lidar clock regularity
    lidar_diffs = diff(Float64.(Dates.datetime2epochms.(win.stare_dt))) ./ 1000.0
    lidar_diffs = lidar_diffs[isfinite.(lidar_diffs) .& (lidar_diffs .> 0)]
    lidar_dt_s     = isempty(lidar_diffs) ? NaN : median(lidar_diffs)
    lidar_dt_std_s = isempty(lidar_diffs) ? NaN : std(lidar_diffs)

    # VN clock regularity inside this chunk window
    ind_vn, vndt_chunk, vn2_chunk = LidarVNSync.vn_subset_20hz(
        win.stare_dt, Vn; pad=Second(0))
    vn_diffs = diff(Float64.(Dates.datetime2epochms.(vndt_chunk))) ./ 1000.0
    vn_diffs = vn_diffs[isfinite.(vn_diffs) .& (vn_diffs .> 0)]
    vn_dt_s     = isempty(vn_diffs) ? NaN : median(vn_diffs)
    vn_dt_std_s = isempty(vn_diffs) ? NaN : std(vn_diffs)
    vn_gap_max_s = isempty(vn_diffs) ? NaN : maximum(vn_diffs)

    # apply stored offset, compute quality metrics
    offset_s = get(offset_by_ic, ic, NaN)
    dt1 = LidarVNSync.dt_seconds(win.stare_dt)
    (!isfinite(dt1) || dt1 <= 0) && (dt1 = LidarVNSync.TIMESTEP)

    corr = NaN; coh = NaN
    if isfinite(offset_s) && offset_s != LidarVNSync.OFFSET_SENTINEL_S
        vn2_shifted = LidarVNSync.shift_signal_linear(win.vn2_xcorr, offset_s; dt=dt1)
        corr = pearson_finite(win.mdv, vn2_shifted)
        coh  = band_coherence_welch_5_20s(win.mdv, vn2_shifted; dt=dt1)
    end

    new_file = (win.ifile != prev_ifile) && (prev_ifile != -1)
    prev_ifile = win.ifile

    # failure diagnosis flags
    flags = String[]
    win.vn_coverage < 0.5 &&
        push!(flags, "LOW_VN_COVERAGE($(round(win.vn_coverage, digits=2)))")
    win.vn2_xcorr_nan_frac > 0.15 &&
        push!(flags, "HIGH_VN_NAN($(round(win.vn2_xcorr_nan_frac, digits=2)))")
    isempty(win.goodlevels) &&
        push!(flags, "NO_MDV_LEVELS")
    count(isfinite, win.mdv) < 30 &&
        push!(flags, "FEW_MDV_FINITE($(count(isfinite, win.mdv)))")
    !isfinite(vn_dt_s) &&
        push!(flags, "NO_VN_DATA")
    isfinite(vn_gap_max_s) && vn_gap_max_s > 2.0 &&
        push!(flags, "VN_GAP($(round(vn_gap_max_s, digits=1))s)")
    isfinite(vn_dt_std_s) && vn_dt_std_s > 0.02 &&
        push!(flags, "VN_CLOCK_JITTER(std=$(round(vn_dt_std_s, digits=3))s)")
    isfinite(lidar_dt_std_s) && lidar_dt_std_s > 0.5 &&
        push!(flags, "LIDAR_CLOCK_JITTER(std=$(round(lidar_dt_std_s, digits=2))s)")
    new_file &&
        push!(flags, "NEW_FILE($(win.ifile))")
    (!isfinite(offset_s) || offset_s == LidarVNSync.OFFSET_SENTINEL_S) &&
        push!(flags, "SENTINEL_OFFSET")
    isfinite(corr) && corr < SYNC_CORR_THR &&
        push!(flags, "LOW_CORR($(round(corr, digits=3)))")
    isfinite(coh)  && coh  < SYNC_COH_THR  &&
        push!(flags, "LOW_COH($(round(coh, digits=3)))")
    isempty(flags) && push!(flags, "ok")
    flag_str = join(flags, "|")

    push!(rows, (;
        ic, ifile=win.ifile, new_file,
        vn_coverage=win.vn_coverage, vn_nan_frac=win.vn2_xcorr_nan_frac,
        ngoodlevels=length(win.goodlevels), mdv_nfinite=count(isfinite, win.mdv),
        lidar_dt_s, lidar_dt_std_s, vn_dt_s, vn_dt_std_s, vn_gap_max_s,
        offset_s, corr, coh, flag=flag_str))

    @printf("%d,%d,%s,%.3f,%.3f,%d,%d,%.4f,%.4f,%.4f,%.4f,%.3f,%.4f,%.4f,%.4f,%s\n",
        ic, win.ifile, new_file,
        win.vn_coverage, win.vn2_xcorr_nan_frac,
        length(win.goodlevels), count(isfinite, win.mdv),
        lidar_dt_s, lidar_dt_std_s,
        isfinite(vn_dt_s)     ? vn_dt_s     : -9999.0,
        isfinite(vn_dt_std_s) ? vn_dt_std_s : -9999.0,
        isfinite(vn_gap_max_s)? vn_gap_max_s: -9999.0,
        isfinite(offset_s)    ? offset_s    : -9999.0,
        isfinite(corr)        ? corr        : -9999.0,
        isfinite(coh)         ? coh         : -9999.0,
        flag_str)
end

# ── save CSV ───────────────────────────────────────────────────────────────

mkpath("epsilon_data")
csv_path = joinpath("epsilon_data", "chunk_failure_scan_$(SCAN_START)_$(SCAN_END).csv")
open(csv_path, "w") do io
    println(io, header)
    for r in rows
        @printf(io, "%d,%d,%s,%.3f,%.3f,%d,%d,%.4f,%.4f,%.4f,%.4f,%.3f,%.4f,%.4f,%.4f,%s\n",
            r.ic, r.ifile, r.new_file,
            r.vn_coverage, r.vn_nan_frac,
            r.ngoodlevels, r.mdv_nfinite,
            r.lidar_dt_s, r.lidar_dt_std_s,
            isfinite(r.vn_dt_s)     ? r.vn_dt_s     : -9999.0,
            isfinite(r.vn_dt_std_s) ? r.vn_dt_std_s : -9999.0,
            isfinite(r.vn_gap_max_s)? r.vn_gap_max_s: -9999.0,
            isfinite(r.offset_s)    ? r.offset_s    : -9999.0,
            isfinite(r.corr)        ? r.corr        : -9999.0,
            isfinite(r.coh)         ? r.coh         : -9999.0,
            r.flag)
    end
end
println("\nSaved: ", csv_path)
