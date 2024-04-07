#
# this file is meant for plotting the global changes in GPP and SIF at different seasons
#
import cartopy.crs as CCRS
import matplotlib.pyplot as PLT
import netCDF4 as NC


# 1. read the data
dset = NC.Dataset("../output/s1_seasonality.nc")
lats = dset.variables["lat"][:]
lons = dset.variables["lon"][:]
dgpp = dset.variables["ΔGPP"][:]
dsif = dset.variables["ΔSIF"][:]
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

# 3. plot the delta GPP and SIF for MAM on ax11 and ax12
ax11.coastlines(linewidth=0.5)
ax11.set_global()
cm11 = ax11.pcolormesh(lons, lats, dgpp[0,:,:], shading="auto", transform=CCRS.PlateCarree(), vmin=-1.5, vmax=1.5, cmap="coolwarm")
fig.colorbar(cm11, ax=ax11, label="gC m$^{-2}$ day$^{-1}$")
ax11.set_title("(a) $\Delta$GPP MAM", fontsize=12, loc="left")

ax12.coastlines(linewidth=0.5)
ax12.set_global()
cm12 = ax12.pcolormesh(lons, lats, dsif[0,:,:], shading="auto", transform=CCRS.PlateCarree(), vmin=-0.15, vmax=0.15, cmap="coolwarm")
fig.colorbar(cm12, ax=ax12, label="mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$")
ax12.set_title("(b) $\Delta$SIF MAM", fontsize=12, loc="left")

# 4. plot the delta GPP and SIF for JJA on ax21 and ax22
ax21.coastlines(linewidth=0.5)
ax21.set_global()
cm21 = ax21.pcolormesh(lons, lats, dgpp[1,:,:], shading="auto", transform=CCRS.PlateCarree(), vmin=-1.5, vmax=1.5, cmap="coolwarm")
fig.colorbar(cm21, ax=ax21, label="gC m$^{-2}$ day$^{-1}$")
ax21.set_title("(c) $\Delta$GPP JJA", fontsize=12, loc="left")

ax22.coastlines(linewidth=0.5)
ax22.set_global()
cm22 = ax22.pcolormesh(lons, lats, dsif[1,:,:], shading="auto", transform=CCRS.PlateCarree(), vmin=-0.15, vmax=0.15, cmap="coolwarm")
fig.colorbar(cm22, ax=ax22, label="mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$")
ax22.set_title("(d) $\Delta$SIF JJA", fontsize=12, loc="left")

# 5. plot the delta GPP and SIF for SON on ax31 and ax32
ax31.coastlines(linewidth=0.5)
ax31.set_global()
cm31 = ax31.pcolormesh(lons, lats, dgpp[2,:,:], shading="auto", transform=CCRS.PlateCarree(), vmin=-1.5, vmax=1.5, cmap="coolwarm")
fig.colorbar(cm31, ax=ax31, label="gC m$^{-2}$ day$^{-1}$")
ax31.set_title("(e) $\Delta$GPP SON", fontsize=12, loc="left")

ax32.coastlines(linewidth=0.5)
ax32.set_global()
cm32 = ax32.pcolormesh(lons, lats, dsif[2,:,:], shading="auto", transform=CCRS.PlateCarree(), vmin=-0.15, vmax=0.15, cmap="coolwarm")
fig.colorbar(cm32, ax=ax32, label="mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$")
ax32.set_title("(f) $\Delta$SIF SON", fontsize=12, loc="left")

# 6. plot the delta GPP and SIF for DJF on ax41 and ax42
ax41.coastlines(linewidth=0.5)
ax41.set_global()
cm41 = ax41.pcolormesh(lons, lats, dgpp[3,:,:], shading="auto", transform=CCRS.PlateCarree(), vmin=-1.5, vmax=1.5, cmap="coolwarm")
fig.colorbar(cm41, ax=ax41, label="gC m$^{-2}$ day$^{-1}$")
ax41.set_title("(g) $\Delta$GPP DJF", fontsize=12, loc="left")

ax42.coastlines(linewidth=0.5)
ax42.set_global()
cm42 = ax42.pcolormesh(lons, lats, dsif[3,:,:], shading="auto", transform=CCRS.PlateCarree(), vmin=-0.15, vmax=0.15, cmap="coolwarm")
fig.colorbar(cm42, ax=ax42, label="mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$")
ax42.set_title("(h) $\Delta$SIF DJF", fontsize=12, loc="left")

# save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/s1_seasonality.png", bbox_inches="tight")
