#
# this script is meant to plot the albedo of the PAR and UV
#
import cartopy.crs as CCRS
import matplotlib.pyplot as PLT
import netCDF4 as NC

# 1. read the data
albedo_ds = NC.Dataset("../output/10_albedo.nc")
lats = albedo_ds.variables["lat"][:]
lons = albedo_ds.variables["lon"][:]
alb_par = albedo_ds.variables["ALB_PAR"][:]
alb_uv = albedo_ds.variables["ALB_UV"][:]
albedo_ds.close()

# 2. create figure 9 with 3 rows
fig = PLT.figure(10, dpi=300, figsize=(6.5,9))
ax1 = fig.add_subplot(3,1,1, projection=CCRS.Robinson())
ax2 = fig.add_subplot(3,1,2, projection=CCRS.Robinson())
ax3 = fig.add_subplot(3,1,3, projection=CCRS.Robinson())
sifig = PLT.figure("10-si", dpi=300, figsize=(6.5,9))
siax1 = sifig.add_subplot(3,1,1, projection=CCRS.Robinson())
siax2 = sifig.add_subplot(3,1,2, projection=CCRS.Robinson())
siax3 = sifig.add_subplot(3,1,3, projection=CCRS.Robinson())

# 3. plot the PAR albedo on ax1
figs = [fig, sifig]
for i in range(2):
    fg = figs[i]
    ax = [ax1,siax1][i]
    ax.coastlines(linewidth=0.5)
    ax.set_global()
    cm = ax.pcolormesh(lons, lats, alb_par, transform=CCRS.PlateCarree(), vmin=0.0, vmax=0.1)
    fg.colorbar(cm, ax=ax)
    ax.set_title("$\\bf{(a)}$ PAR albedo", fontsize=12, loc="left")

# 4. plot the UV albedo on ax2
for i in range(2):
    fg = figs[i]
    ax = [ax2,siax2][i]
    ax.coastlines(linewidth=0.5)
    ax.set_global()
    cm = ax.pcolormesh(lons, lats, alb_uv, transform=CCRS.PlateCarree(), vmin=0.0, vmax=0.1)
    fg.colorbar(cm, ax=ax)
    ax.set_title("$\\bf{(b)}$ UV albedo", fontsize=12, loc="left")

# 5. plot the difference between PAR and UV albedo on ax3
vmaxs = [0.1, 0.02]
for i in range(2):
    fg = figs[i]
    ax = [ax3,siax3][i]
    ax.coastlines(linewidth=0.5)
    ax.set_global()
    cm = ax.pcolormesh(lons, lats, alb_par - alb_uv, transform=CCRS.PlateCarree(), vmin=0.0, vmax=vmaxs[i])
    fg.colorbar(cm, ax=ax)
    ax.set_title("$\\bf{(c)}$ PAR $-$ UV albedo", fontsize=12, loc="left")

# save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/10_albedo.png", bbox_inches="tight")
sifig.set_tight_layout(True)
sifig.savefig("../figure/s8_albedo.png", bbox_inches="tight")
