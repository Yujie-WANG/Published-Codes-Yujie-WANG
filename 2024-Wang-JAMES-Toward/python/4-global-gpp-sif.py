#
# this file is meant for plotting the global changes in GPP and SIF at different scenarios
#
import cartopy.crs as CCRS
import matplotlib.pyplot as PLT
import netCDF4 as NC


# 1. read the data
dset = NC.Dataset("../output/4_gpp_sif.nc")
lats = dset.variables["lat"][:]
lons = dset.variables["lon"][:]
gpp_hs_kdtd = dset.variables["GPP_HS_KDTD"][:]
gpp_hs_kdst = dset.variables["GPP_HS_KDST"][:]
gpp_bb_kdtd = dset.variables["GPP_BB_KDTD"][:]
gpp_bb_kdst = dset.variables["GPP_BB_KDST"][:]
sif_hs_kdtd = dset.variables["SIF_HS_KDTD"][:]
sif_hs_kdst = dset.variables["SIF_HS_KDST"][:]
sif_bb_kdtd = dset.variables["SIF_BB_KDTD"][:]
sif_bb_kdst = dset.variables["SIF_BB_KDST"][:]
dset.close()

# 2. create figure 4 with 4 rows and 2 columns
fig = PLT.figure(4, dpi=300, figsize=(12,11))
ax11 = fig.add_subplot(4,2,1, projection=CCRS.Robinson())
ax12 = fig.add_subplot(4,2,2, projection=CCRS.Robinson())
ax21 = fig.add_subplot(4,2,3, projection=CCRS.Robinson())
ax22 = fig.add_subplot(4,2,4, projection=CCRS.Robinson())
ax31 = fig.add_subplot(4,2,5, projection=CCRS.Robinson())
ax32 = fig.add_subplot(4,2,6, projection=CCRS.Robinson())
ax41 = fig.add_subplot(4,2,7, projection=CCRS.Robinson())
ax42 = fig.add_subplot(4,2,8, projection=CCRS.Robinson())

# 3. plot the delta GPP and SIF for HS_KDTD and HS_KDST on ax11 and ax12
ax11.coastlines(linewidth=0.5)
ax11.set_global()
cm11 = ax11.pcolormesh(lons, lats, gpp_hs_kdst - gpp_hs_kdtd, shading="auto", transform=CCRS.PlateCarree(), vmin=-0.15, vmax=0.15, cmap="coolwarm")
fig.colorbar(cm11, ax=ax11, label="gC m$^{-2}$ day$^{-1}$")
ax11.set_title("(a) $\Delta$GPP HS_KDST $-$ HS_KDTD", fontsize=12, loc="left")

ax12.coastlines(linewidth=0.5)
ax12.set_global()
cm12 = ax12.pcolormesh(lons, lats, sif_hs_kdst - sif_hs_kdtd, shading="auto", transform=CCRS.PlateCarree(), vmin=-0.04, vmax=0.04, cmap="coolwarm")
fig.colorbar(cm12, ax=ax12, label="mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$")
ax12.set_title("(b) $\Delta$SIF HS_KDST $-$ HS_KDTD", fontsize=12, loc="left")

# 4. plot the delta GPP and SIF for HS_NVER and BB_KDTD on ax21 and ax22
ax21.coastlines(linewidth=0.5)
ax21.set_global()
cm21 = ax21.pcolormesh(lons, lats, gpp_bb_kdtd - gpp_hs_kdtd, shading="auto", transform=CCRS.PlateCarree(), vmin=-1.2, vmax=1.2, cmap="coolwarm")
fig.colorbar(cm21, ax=ax21, label="gC m$^{-2}$ day$^{-1}$")
ax21.set_title("(c) $\Delta$GPP BB_KDTD $-$ HS_KDTD", fontsize=12, loc="left")

ax22.coastlines(linewidth=0.5)
ax22.set_global()
cm22 = ax22.pcolormesh(lons, lats, sif_bb_kdtd - sif_hs_kdtd, shading="auto", transform=CCRS.PlateCarree(), vmin=-0.1, vmax=0.1, cmap="coolwarm")
fig.colorbar(cm22, ax=ax22, label="mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$")
ax22.set_title("(d) $\Delta$SIF BB_KDTD $-$ HS_KDTD", fontsize=12, loc="left")

# 5. plot the delta GPP and SIF for HS_KDTD and BB_KDST on ax31 and ax32
ax31.coastlines(linewidth=0.5)
ax31.set_global()
cm31 = ax31.pcolormesh(lons, lats, gpp_bb_kdst - gpp_hs_kdtd, shading="auto", transform=CCRS.PlateCarree(), vmin=-1.2, vmax=1.2, cmap="coolwarm")
fig.colorbar(cm31, ax=ax31, label="gC m$^{-2}$ day$^{-1}$")
ax31.set_title("(e) $\Delta$GPP BB_KDST $-$ HS_KDTD", fontsize=12, loc="left")

ax32.coastlines(linewidth=0.5)
ax32.set_global()
cm32 = ax32.pcolormesh(lons, lats, sif_bb_kdst - sif_hs_kdtd, shading="auto", transform=CCRS.PlateCarree(), vmin=-0.12, vmax=0.12, cmap="coolwarm")
fig.colorbar(cm32, ax=ax32, label="mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$")
ax32.set_title("(f) $\Delta$SIF BB_KDST $-$ HS_KDTD", fontsize=12, loc="left")

# 6. plot the relative change in GPP and SIF for HS_KDTD and BB_KDST on ax41 and ax42
ax41.coastlines(linewidth=0.5)
ax41.set_global()
cm41 = ax41.pcolormesh(lons, lats, gpp_bb_kdst / gpp_hs_kdtd, shading="auto", transform=CCRS.PlateCarree(), vmin=0.8, vmax=1.2, cmap="coolwarm")
fig.colorbar(cm41, ax=ax41)
ax41.set_title("(g) GPP BB_KDST / HS_KDTD", fontsize=12, loc="left")

ax42.coastlines(linewidth=0.5)
ax42.set_global()
cm42 = ax42.pcolormesh(lons, lats, sif_bb_kdst / sif_hs_kdtd, shading="auto", transform=CCRS.PlateCarree(), vmin=0.7, vmax=1.3, cmap="coolwarm")
fig.colorbar(cm42, ax=ax42)
ax42.set_title("(h) SIF BB_KDST / HS_KDTD", fontsize=12, loc="left")

# save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/4_gpp_sif.png", bbox_inches="tight")
