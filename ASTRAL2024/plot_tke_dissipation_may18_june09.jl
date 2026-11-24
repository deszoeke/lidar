using Pkg
Pkg.activate(@__DIR__)

using Dates
using JLD2
using PythonPlot
using Statistics

include(joinpath(@__DIR__, "compute_hmix_epsilon.jl"))

data = load(joinpath(@__DIR__, "epsilon_data", "epsi_stare_vn_chunks_r20260907_1008.jld"))
epsi = data["epsi"]
height = data["height"]
start_dt = data["dtime_st"]

plot_epsi = similar(epsi, Float64)
for i in eachindex(epsi)
    e = epsi[i]
    plot_epsi[i] = (ismissing(e) || !isfinite(e) || e <= 0) ? NaN : e
end

# Temporary visualization only: the old MDV file's post-gap times are six days early.
mdv = load(joinpath(@__DIR__, "epsilon_data", "epsi_stare_chunks_vn0_r20260412_1306.jld"))
isequal(mdv["height"], height) || error("MDV and current height grids differ")
issorted(start_dt) || error("Current chunk times are not sorted")
matched_rows = Set{Int}()
unmatched_rows = 0
for (i, old_time) in enumerate(mdv["dtime_st"])
    shifted_time = old_time + Day(6)
    DateTime(2024, 5, 19) <= shifted_time < DateTime(2024, 6, 1) || continue
    j = searchsortedfirst(start_dt, shifted_time)
    candidates = filter(k -> 1 <= k <= length(start_dt), [j - 1, j])
    k = candidates[argmin(abs.(Dates.value.(start_dt[candidates] .- shifted_time)))]
    if abs(Dates.value(start_dt[k] - shifted_time)) > 10
        global unmatched_rows += 1
        continue
    end
    k in matched_rows && error("Multiple MDV chunks match current chunk $k")
    push!(matched_rows, k)
    for iz in axes(plot_epsi, 2)
        e = mdv["epsi"][i, iz]
        if !isfinite(plot_epsi[k, iz]) && !ismissing(e) && isfinite(e) && e > 0
            plot_epsi[k, iz] = e
        end
    end
end
unmatched_rows == 0 || @warn "Provisional MDV chunks did not match current timestamps" unmatched_rows
println("Matched $(length(matched_rows)) provisional MDV chunks")

# Keep long gaps blank rather than stretching an adjacent pcolormesh cell across them.
mesh_time = DateTime[]
mesh_rows = Vector{Vector{Float64}}()
for i in eachindex(start_dt)
    push!(mesh_time, start_dt[i])
    push!(mesh_rows, plot_epsi[i, 4:end])
    if i < length(start_dt) && start_dt[i + 1] - start_dt[i] > Minute(30)
        push!(mesh_time, start_dt[i] + Minute(10), start_dt[i + 1] - Minute(10))
        push!(mesh_rows, fill(NaN, length(height) - 3), fill(NaN, length(height) - 3))
    end
end

PythonPlot.matplotlib.rcParams.update(Dict(
    "font.family" => "sans-serif",
    "font.sans-serif" => ["Arial", "Helvetica", "OpenSans","DejaVu Sans"],
    "font.size" => 14,
    "axes.titlesize" => 17,
    "axes.labelsize" => 16,
    "xtick.labelsize" => 13,
    "ytick.labelsize" => 13,
    "legend.fontsize" => 13,
))
fig = figure(figsize=(13, 5.5))
norm = PythonPlot.matplotlib.colors.LogNorm(vmin=1e-6, vmax=1e-3)
pc = pcolormesh(mesh_time, height[4:end] ./ 1e3, reduce(hcat, mesh_rows),
            cmap="RdYlBu_r", norm=norm)
gca().set_facecolor(0.6 .+ [0, 0, 0])
gca().xaxis.set_major_formatter(PythonPlot.matplotlib.dates.DateFormatter("%m%d"))
xlim([DateTime(2024, 5, 18), DateTime(2024, 6, 10)])
ylim([0, 1.5])
ylabel("height (km)")
xlabel("date 2024")
title("EKAMSAT TKE dissipation rate: provisional MDV before June 1")
axvline(DateTime(2024, 6, 1), color="k", linestyle="--", linewidth=1)
cb = colorbar(pc, pad=0.12)
cb.set_label(L"\epsilon\;(\mathrm{m}^2\ \mathrm{s}^{-3})")

_, hmix_height_m, _ = compute_hmix_for_matrix(plot_epsi[:, 4:end], height[4:end])
window = (start_dt .>= DateTime(2024, 5, 18)) .& (start_dt .< DateTime(2024, 6, 10))
plot_time = start_dt[window]
hmix_raw = [ismissing(h) ? NaN : h / 1e3 for h in hmix_height_m[window]]
hmix_plot = fill(NaN, length(hmix_raw))
for i in eachindex(hmix_raw)
    isfinite(hmix_raw[i]) || continue
    lo = hi = i
    while lo > 1 && isfinite(hmix_raw[lo - 1]) &&
          plot_time[i] - plot_time[lo - 1] <= Hour(1) &&
          plot_time[lo] - plot_time[lo - 1] <= Minute(30)
        lo -= 1
    end
    while hi < length(hmix_raw) && isfinite(hmix_raw[hi + 1]) &&
          plot_time[hi + 1] - plot_time[i] <= Hour(1) &&
          plot_time[hi + 1] - plot_time[hi] <= Minute(30)
        hi += 1
    end
    hmix_plot[i] = mean(@view hmix_raw[lo:hi])
end
ax_hmix = twinx()
plot(plot_time, hmix_plot, color="k", linewidth=1.0, alpha=0.85,
     label="h_mix (2 h mean)")
ax_hmix.set_ylim([0, 1.5])
ax_hmix.set_ylabel("mixed-layer height (km)", color="k")
ax_hmix.tick_params(axis="y", labelcolor="k")
# ax_hmix.legend(loc="upper right")
tight_layout()

outfile = joinpath(@__DIR__, "epsilon_data", "tke_dissipation_may18_june09_provisional_mdv.png")
savefig(outfile, dpi=300)
println(outfile)
