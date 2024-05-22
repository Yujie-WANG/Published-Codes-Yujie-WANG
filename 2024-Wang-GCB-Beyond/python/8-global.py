#
# this script is meant to plot the global sensitivity analysis of the model
#
import cartopy.crs as CCRS
import matplotlib.pyplot as PLT
import netCDF4 as NC


# 1. read the data
dset = NC.Dataset("../output/8_global.nc")
lat = dset.variables["lat"][:]
lon = dset.variables["lon"][:]
dgpp = dset.variables["ΔGPP"][:]
sgpp = dset.variables["SΔGPP"][:]
rgpp = dset.variables["RGPP"][:]
mlai = dset.variables["mLAI"][:]
dset.close()

# 2. create figure 7 with 1 row and 2 columns (the with ratio is 2:1)
fig = PLT.figure(7, dpi=300, figsize=(10.5,6.5))
fig.subplots(nrows=2, ncols=2, width_ratios=[2.4,1], subplot_kw={"projection": CCRS.Robinson()})
axs = fig.axes
axs[1].remove()
axs[1] = fig.add_subplot(222)
axs[3].remove()
axs[3] = fig.add_subplot(224)

# 3. plot the data
axs[0].coastlines(linewidth=0.5)
axs[0].set_global()
cm = axs[0].pcolormesh(lon, lat, dgpp, transform=CCRS.PlateCarree(), vmin=0, vmax=0.6)
cb = PLT.colorbar(cm, ax=axs[0], fraction=0.05, pad=0.05)
axs[0].set_title("$\\bf{(a)}$", fontsize=12, loc="left")

mask = dgpp != 0
cm = axs[1].hexbin(mlai[mask], dgpp[mask], gridsize=50, bins="log", cmap="Greys", mincnt=1)
cb = PLT.colorbar(cm, ax=axs[1], fraction=0.05, pad=0.05)
axs[1].text(0.05, 0.6, "R$^2$ = 0.873", fontsize=12)
axs[1].set_ylabel("$\Delta$GPP (μmol m$^{-2}$ s$^{-1}$)", fontsize=12)
axs[1].set_title("$\\bf{(b)}$", fontsize=12, loc="left")

axs[2].coastlines(linewidth=0.5)
axs[2].set_global()
cm = axs[2].pcolormesh(lon, lat, rgpp, transform=CCRS.PlateCarree(), vmin=0, vmax=6)
cb = PLT.colorbar(cm, ax=axs[2], fraction=0.05, pad=0.05)
axs[2].set_title("$\\bf{(c)}$", fontsize=12, loc="left")

mask = rgpp != 0
cm = axs[3].hexbin(mlai[mask], rgpp[mask], gridsize=50, bins="log", cmap="Greys", mincnt=1)
cb = PLT.colorbar(cm, ax=axs[3], fraction=0.05, pad=0.05)
axs[3].text(0.05, 6, "R$^2$ = 0.867", fontsize=12)
axs[3].set_xlabel("Annual mean LAI (-)", fontsize=12)
axs[3].set_ylabel("Relative $\Delta$GPP (%)", fontsize=12)
axs[3].set_title("$\\bf{(d)}$", fontsize=12, loc="left")

# 4. save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/8_global.png", bbox_inches="tight")

# 5. create a si figure s5 with two rows and two columns
sifig = PLT.figure("s5", dpi=300, figsize=(12,6))
sifig.subplots(nrows=2, ncols=2, subplot_kw={"projection": CCRS.Robinson()})
siaxs = sifig.axes

siaxs[0].coastlines(linewidth=0.5)
siaxs[0].set_global()
cb = siaxs[0].pcolormesh(lon, lat, sgpp[0,:,:], transform=CCRS.PlateCarree(), vmin=0, vmax=0.6)
PLT.colorbar(cb, ax=siaxs[0], fraction=0.05, pad=0.05)
siaxs[0].set_title("(a) MAM", fontsize=12, loc="left")

siaxs[1].coastlines(linewidth=0.5)
siaxs[1].set_global()
cb = siaxs[1].pcolormesh(lon, lat, sgpp[1,:,:], transform=CCRS.PlateCarree(), vmin=0, vmax=0.6)
PLT.colorbar(cb, ax=siaxs[1], fraction=0.05, pad=0.05)
siaxs[1].set_title("(b) JJA", fontsize=12, loc="left")

siaxs[2].coastlines(linewidth=0.5)
siaxs[2].set_global()
cb = siaxs[2].pcolormesh(lon, lat, sgpp[2,:,:], transform=CCRS.PlateCarree(), vmin=0, vmax=0.6)
PLT.colorbar(cb, ax=siaxs[2], fraction=0.05, pad=0.05)
siaxs[2].set_title("(c) SON", fontsize=12, loc="left")

siaxs[3].coastlines(linewidth=0.5)
siaxs[3].set_global()
cb = siaxs[3].pcolormesh(lon, lat, sgpp[3,:,:], transform=CCRS.PlateCarree(), vmin=0, vmax=0.6)
PLT.colorbar(cb, ax=siaxs[3], fraction=0.05, pad=0.05)
siaxs[3].set_title("(d) DJF", fontsize=12, loc="left")

sifig.set_tight_layout(True)
sifig.savefig("../figure/s7_seasonality.png", bbox_inches="tight")
