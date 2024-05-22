#
# this script is meant to plot the model testing results
#
import matplotlib.pyplot as PLT
import numpy as NP
import pandas as PD


# 1. read in the data
fitting = PD.read_csv("../output/6_fitting.csv")

# 2. create a figure to plot on
fig = PLT.figure(4, dpi=300, figsize=(6,3.3))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# 3. plot the data comparison on ax1 and ax2 (points)
mask_fr = fitting.LIGHT == "FR"
mask_par = fitting.LIGHT == "PAR"
mask_combo = fitting.LIGHT == "PARFR"
ax1.plot(fitting.GPP[mask_fr], fitting.GPP_FIT[mask_fr], "rx", label="FR")
ax1.plot(fitting.GPP[mask_par], fitting.GPP_FIT[mask_par], "go", mfc="none", label="PAR")
ax1.plot(fitting.GPP[mask_combo], fitting.GPP_FIT[mask_combo], "b+", label="PAR+FR")
ax2.plot(fitting.GPP[mask_fr], fitting.GPP_FIT_CYTO[mask_fr], "rx", label="FR")
ax2.plot(fitting.GPP[mask_par], fitting.GPP_FIT_CYTO[mask_par], "go", mfc="none", label="PAR")
ax2.plot(fitting.GPP[mask_combo], fitting.GPP_FIT_CYTO[mask_combo], "b+", label="PAR+FR")
for ax in [ax1, ax2]:
    ax.plot([-1, 36], [-1, 36], "k--")
    ax.set_xlim(-1, 36)
    ax.set_ylim(-1, 36)
    ax.set_xlabel("obs GPP (µmol m$^{-2}$ s$^{-1}$)", fontsize=12)
ax1.legend(loc="best")
ax1.set_ylabel("mod GPP (µmol m$^{-2}$ s$^{-1}$)", fontsize=12)
ax1.set_title("$\\bf{(a)}$ Farquhar et al. model", fontsize=12, loc="left")
ax2.set_title("$\\bf{(b)}$ Johnson & Berry model", fontsize=12, loc="left")

# 4. plot the RMSE on ax1
def rmse(predictions, targets):
    return NP.sqrt(((predictions - targets) ** 2).mean())

ax1.text(35, 0, "RMSE = %.2f" % rmse(fitting.GPP[mask_fr], fitting.GPP_FIT[mask_fr]), color="r", ha="right", va="bottom")
ax1.text(35, 4, "RMSE = %.2f" % rmse(fitting.GPP[mask_combo], fitting.GPP_FIT[mask_combo]), color="b", ha="right", va="bottom")
ax1.text(35, 8, "RMSE = %.2f" % rmse(fitting.GPP[mask_par], fitting.GPP_FIT[mask_par]), color="g", ha="right", va="bottom")
ax2.text(35, 0, "RMSE = %.2f" % rmse(fitting.GPP[mask_fr], fitting.GPP_FIT_CYTO[mask_fr]), color="r", ha="right", va="bottom")
ax2.text(35, 4, "RMSE = %.2f" % rmse(fitting.GPP[mask_combo], fitting.GPP_FIT_CYTO[mask_combo]), color="b", ha="right", va="bottom")
ax2.text(35, 8, "RMSE = %.2f" % rmse(fitting.GPP[mask_par], fitting.GPP_FIT_CYTO[mask_par]), color="g", ha="right", va="bottom")

# 5. save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/6_model_testing.pdf", bbox_inches="tight")

# 6. plot the results per TREE in a 4*8 figure
sifig = PLT.figure("S3", dpi=300, figsize=(8,12))
sifig_cyto = PLT.figure("S4", dpi=300, figsize=(8,12))
axs = []
axs_cyto = []
for i in range(8):
    for j in range(4):
        axs.append(sifig.add_subplot(8,4,i*4+j+1))
        axs_cyto.append(sifig_cyto.add_subplot(8,4,i*4+j+1))

for i in range(32):
    subdf = fitting[fitting.TREE == i+1]
    mask_fr = subdf.LIGHT == "FR"
    mask_par = subdf.LIGHT == "PAR"
    mask_par_fr = subdf.LIGHT == "PARFR"
    axs[i].plot(subdf.TPAR[mask_fr], subdf.GPP[mask_fr], "r+", label="FR")
    axs[i].plot(subdf.TPAR[mask_fr], subdf.GPP_FIT[mask_fr], "ro", mfc="none")
    axs[i].plot(subdf.TPAR[mask_par], subdf.GPP[mask_par], "g+", label="PAR")
    axs[i].plot(subdf.TPAR[mask_par], subdf.GPP_FIT[mask_par], "go", mfc="none")
    axs[i].plot(subdf.TPAR[mask_par_fr], subdf.GPP[mask_par_fr], "b+", label="PAR+FR")
    axs[i].plot(subdf.TPAR[mask_par_fr], subdf.GPP_FIT[mask_par_fr], "bo", mfc="none")
    axs_cyto[i].plot(subdf.TPAR[mask_fr], subdf.GPP[mask_fr], "r+", label="FR")
    axs_cyto[i].plot(subdf.TPAR[mask_fr], subdf.GPP_FIT_CYTO[mask_fr], "ro", mfc="none")
    axs_cyto[i].plot(subdf.TPAR[mask_par], subdf.GPP[mask_par], "g+", label="PAR")
    axs_cyto[i].plot(subdf.TPAR[mask_par], subdf.GPP_FIT_CYTO[mask_par], "go", mfc="none")
    axs_cyto[i].plot(subdf.TPAR[mask_par_fr], subdf.GPP[mask_par_fr], "b+", label="PAR+FR")
    axs_cyto[i].plot(subdf.TPAR[mask_par_fr], subdf.GPP_FIT_CYTO[mask_par_fr], "bo", mfc="none")

for i in [0,4,8,12,16,20,24,28]:
    axs[i].set_ylabel("GPP\n(µmol m$^{-2}$ s$^{-1}$)", fontsize=10)
    axs_cyto[i].set_ylabel("GPP\n(µmol m$^{-2}$ s$^{-1}$)", fontsize=10)

for i in [28,29,30,31]:
    axs[i].set_xlabel("PAR (µmol m$^{-2}$ s$^{-1}$)", fontsize=10)
    axs_cyto[i].set_xlabel("PAR (µmol m$^{-2}$ s$^{-1}$)", fontsize=10)

sifig.set_tight_layout(True)
sifig.savefig("../figure/s4_model_testing.pdf", bbox_inches="tight")
sifig_cyto.set_tight_layout(True)
sifig_cyto.savefig("../figure/s5_model_testing_cyto.pdf", bbox_inches="tight")
