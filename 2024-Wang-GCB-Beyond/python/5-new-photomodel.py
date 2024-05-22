#
# this script is meant to plot the model output to illustrate that the new model is working as expected
#
import pandas as PD
import matplotlib.pyplot as PLT


# 1. read the data
sensitivity = PD.read_csv("../output/5_fr_diff.csv")
profiles = PD.read_csv("../output/5_fr_diff_profiles.csv")

# 2. create a figure to plot on
fig = PLT.figure(3, dpi=300, figsize=(6,6))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(234)
ax3 = fig.add_subplot(235)
ax4 = fig.add_subplot(236)
sifig = PLT.figure("S2", dpi=300, figsize=(6,6))
siax1 = sifig.add_subplot(211)
siax2 = sifig.add_subplot(234)
siax3 = sifig.add_subplot(235)
siax4 = sifig.add_subplot(236)

# 3. plot the A~PAR relationship on ax1
ax1.plot(sensitivity.Ratio, sensitivity.RPAR_UV, "k:", label="UV+PAR")
ax1.plot(sensitivity.Ratio, sensitivity.RPAR_UV_FR, "k-", label="UV+PAR+FR")
ax1.plot(sensitivity.Ratio, sensitivity.FR, "r:", label="FR")
ax1.plot(sensitivity.Ratio, sensitivity["2%RPAR_FR"], "r-", label="FR+2%UV+2%PAR")
siax1.plot(sensitivity.Ratio, sensitivity.RPAR_UV_CYTO, "k:", label="UV+PAR")
siax1.plot(sensitivity.Ratio, sensitivity.RPAR_UV_FR_CYTO, "k-", label="UV+PAR+FR")
siax1.plot(sensitivity.Ratio, sensitivity.FR_CYTO, "r:", label="FR")
siax1.plot(sensitivity.Ratio, sensitivity["2%RPAR_FR_CYTO"], "r-", label="FR+2%UV+2%PAR")
for ax in [ax1, siax1]:
    ax.legend(loc="center right")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.1, 27)
    ax.set_xlabel("Ratio of Natural Radiation (-)", fontsize=12)
    ax.set_ylabel("GPP (µmol m$^{-2}$ s$^{-1}$)", fontsize=12)
    ax.set_title("$\\bf{(a)}$", fontsize=12, loc="left")

# 4. plot the profiles on ax2, ax3, ax4
ax2.plot(profiles.UV, list(reversed(range(1,13))), "-", color="darkviolet", label="UV")
ax2.plot(profiles.PAR, list(reversed(range(1,13))), "-", color="forestgreen", label="PAR")
ax2.plot(profiles.FR, list(reversed(range(1,13))), "-", color="red", label="FR")
siax2.plot(profiles.UV_CYTO, list(reversed(range(1,13))), "-", color="darkviolet", label="UV")
siax2.plot(profiles.PAR_CYTO, list(reversed(range(1,13))), "-", color="forestgreen", label="PAR")
siax2.plot(profiles.FR_CYTO, list(reversed(range(1,13))), "-", color="red", label="FR")
for ax in [ax2, siax2]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, 450)
    ax.set_ylim(0, 12)
    ax.legend(loc="best")
    ax.set_xlabel("Rad (W m$^{-2}$)", fontsize=12)
    ax.set_ylabel("Canopy Layer from Bottom", fontsize=12)
    ax.set_title("$\\bf{(b)}$", fontsize=12, loc="left")

ax3.plot(profiles.AG_750, list(reversed(range(1,13))), "k-", label="750 nm")
ax3.plot(profiles.AG_700, list(reversed(range(1,13))), "k:", label="700 nm")
siax3.plot(profiles.AG_750_CYTO, list(reversed(range(1,13))), "k-", label="750 nm")
siax3.plot(profiles.AG_700_CYTO, list(reversed(range(1,13))), "k:", label="700 nm")
for ax in [ax3, siax3]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(3, 14)
    ax.set_ylim(0, 12)
    ax.set_xlabel("$A_g$ (μmol m$^{-2}$ s$^{-1}$)", fontsize=12)
    ax.set_title("$\\bf{(c)}$", fontsize=12, loc="left")

ax4.plot(profiles.GS_750, list(reversed(range(1,13))), "k-", label="750 nm")
ax4.plot(profiles.GS_700, list(reversed(range(1,13))), "k:", label="700 nm")
siax4.plot(profiles.GS_750_CYTO, list(reversed(range(1,13))), "k-", label="750 nm")
siax4.plot(profiles.GS_700_CYTO, list(reversed(range(1,13))), "k:", label="700 nm")
for ax in [ax4, siax4]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, 12)
    ax.legend(loc="best")
    ax.set_xlabel("$g_s$ (mol m$^{-2}$ s$^{-1}$)", fontsize=12)
    ax.set_title("$\\bf{(d)}$", fontsize=12, loc="left")

# save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/5_fr_diff.pdf", bbox_inches="tight")
sifig.set_tight_layout(True)
sifig.savefig("../figure/s3_fr_diff.pdf", bbox_inches="tight")
