cd("/Users/deszoeks/Projects/lidar/ASTRAL2024/")

using Pkg; Pkg.activate(".")
using NCDatasets
using Statistics
using PythonPlot

rc("font", family="sans-serif")
rc("font"; Symbol("sans-serif") => "Arial")
np = PythonPlot.pyimport("numpy")

uvdir = "./data/netcdf/v1"
UV = NCDataset(joinpath(uvdir, "ekamsat_lidar_uv_20240429-20240613_v1.nc"))

# for plotting masked arrays
snr = UV["snr"][:, :]
thr = 1.03
ut_masked = np.ma.masked_array(UV["ut"][:, :], mask=(snr .< thr))
vt_masked = np.ma.masked_array(UV["vt"][:, :], mask=(snr .< thr))

fig, axs = subplots(2,1, sharex=true, figsize=(10,6))
ax = axs[0]
mesh = ax.pcolormesh(UV["time"][:], UV["z_range"][:]/1e3, ut_masked, 
    shading="nearest", cmap="RdBu_r", vmin= -10, vmax=10)
colorbar(mesh, ax=ax, label="U (m/s)")
ax.set_facecolor("lightgray")
ax.set_ylabel("height (km)")

ax = axs[1]
mesh = ax.pcolormesh(UV["time"][:], UV["z_range"][:]/1e3, vt_masked, 
    shading="nearest", cmap="RdBu_r", vmin= -10, vmax=10)
colorbar(mesh, ax=ax, label="V (m/s)")
ax.set_facecolor("lightgray")
ax.set_ylabel("height (km)")
ax.set_xlabel("time")

tight_layout()
display(fig)
