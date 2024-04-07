#
# this script is meant to analyze how the global GPP and SIF biases are related to PFT
#
import matplotlib.pyplot as PLT
import netCDF4 as NC
import numpy as NP


# 1. read the data
dset_gpp = NC.Dataset("../output/4_gpp_sif.nc")
dset_pft = NC.Dataset("../output/5_pft.nc")
delta_gpp = dset_gpp.variables["GPP_BB_KDST"][:] - dset_gpp.variables["GPP_HS_KDTD"][:]
delta_sif = dset_gpp.variables["SIF_BB_KDST"][:] - dset_gpp.variables["SIF_HS_KDTD"][:]
pft = dset_pft.variables["PFT"][:]
dset_gpp.close()
dset_pft.close()

# 2. compute the mean and std of the delta GPP and SIF for each PFT
mean_gpps = []
std_gpps = []
mean_sifs = []
std_sifs = []
for i in range(2,13):
    mask = (pft==i) & ((delta_gpp > 0) | (delta_gpp < 0))
    mean_gpp = NP.mean(delta_gpp[mask])
    std_gpp = NP.std(delta_gpp[mask])
    mean_sif = NP.mean(delta_sif[mask])
    std_sif = NP.std(delta_sif[mask])
    mean_gpps.append(mean_gpp)
    std_gpps.append(std_gpp)
    mean_sifs.append(mean_sif)
    std_sifs.append(std_sif)

# 3. create figure 5 with 1 row and 2 columns to plot
pfts = ["NET-P", "NET-B", "NDT-B", "BET-T", "BET-P", "BDT-T", "BDT-P", "BDT-B", "BES-P", "BDS-P", "BDS-B"]
fig = PLT.figure(5, dpi=300, figsize=(7.5,3.7))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.bar(range(2,13), mean_gpps, yerr=std_gpps, color="gray")
ax1.set_xticks(range(2,13))
ax1.set_xticklabels(pfts, rotation=90)
ax1.set_ylabel("$\Delta$GPP (gC m$^{-2}$ day$^{-1}$)", fontsize=12)
ax1.set_title("(a)", fontsize=12, loc="left")

ax2.bar(range(2,13), mean_sifs, yerr=std_sifs, color="gray")
ax2.set_xticks(range(2,13))
ax2.set_xticklabels(pfts, rotation=90)
ax2.set_ylabel("$\Delta$SIF (mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$)", fontsize=12)
ax2.set_title("(b)", fontsize=12, loc="left")

# 4. save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/5_pft.pdf", bbox_inches="tight")
