#
# this file is meant to plot the global bias in f_par * f_ppar
#
import cartopy.crs as CCRS
import matplotlib.pyplot as PLT
import netCDF4 as NC


# 1. read the data
bias_ds = NC.Dataset("../output/3_bias.nc")
lats = bias_ds.variables["lat"][:]
lons = bias_ds.variables["lon"][:]
ff_toc = bias_ds.variables["FF_TOC"][:]
ff_boc = bias_ds.variables["FF_BOC"][:]
bias_ds.close()

# 2. create figure 3 with 2 rows and 2 columns
fig = PLT.figure(3, dpi=300, figsize=(12,5.5))
ax1 = fig.add_subplot(2,2,1, projection=CCRS.Robinson())
ax2 = fig.add_subplot(2,2,2, projection=CCRS.Robinson())
ax3 = fig.add_subplot(2,2,3, projection=CCRS.Robinson())
ax4 = fig.add_subplot(2,2,4, projection=CCRS.Robinson())

# 3. plot the TOC fapar * fppar on ax1
ax1.coastlines(linewidth=0.5)
ax1.set_global()
cm1 = ax1.pcolormesh(lons, lats, ff_toc, shading="auto", transform=CCRS.PlateCarree(), vmin=0.0, vmax=0.7)
fig.colorbar(cm1, ax=ax1)
ax1.set_title("(a) TOC $f_\mathrm{APAR} \cdot f_\mathrm{PPAR}$", fontsize=12, loc="left")

# 4. plot the BOC fapar * fppar on ax2
ax2.coastlines(linewidth=0.5)
ax2.set_global()
cm2 = ax2.pcolormesh(lons, lats, ff_boc, shading="auto", transform=CCRS.PlateCarree(), vmin=0.0, vmax=0.7)
fig.colorbar(cm2, ax=ax2)
ax2.set_title("(b) BOC $f_\mathrm{APAR} \cdot f_\mathrm{PPAR}$", fontsize=12, loc="left")

# 5. plot the bias on ax3
ax3.coastlines(linewidth=0.5)
ax3.set_global()
cm3 = ax3.pcolormesh(lons, lats, 0.86 - ff_toc, cmap="coolwarm", shading="auto", transform=CCRS.PlateCarree(), vmin=-0.5, vmax=0.5)
fig.colorbar(cm3, ax=ax3)
ax3.set_title("(c) CLM $-$ CliMA TOC $f_\mathrm{APAR} \cdot f_\mathrm{PPAR}$", fontsize=12, loc="left")

# 6. plot the bias on ax4
ax4.coastlines(linewidth=0.5)
ax4.set_global()
cm4 = ax4.pcolormesh(lons, lats, ff_toc - ff_boc, cmap="coolwarm", shading="auto", transform=CCRS.PlateCarree(), vmin=-0.12, vmax=0.12)
fig.colorbar(cm4, ax=ax4)
ax4.set_title("(d) CliMA TOC $-$ BOC $f_\mathrm{APAR} \cdot f_\mathrm{PPAR}$", fontsize=12, loc="left")

# save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/3_bias.png", bbox_inches="tight")
