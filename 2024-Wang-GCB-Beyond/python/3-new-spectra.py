#
# This script is meant to plot the model output to illustrate that the new model is working as expected
#
import pandas as PD
import matplotlib.pyplot as PLT

# 1. read the data
uvdiff = PD.read_csv("../output/3_uv_diff.csv")

# 2. create a figure to plot on
fig = PLT.figure(2, dpi=300, figsize=(6,3))
ax2 = fig.add_subplot(131)
ax3 = fig.add_subplot(132)
ax4 = fig.add_subplot(133)
sifig = PLT.figure("S1", dpi=300, figsize=(6,3))
siax2 = sifig.add_subplot(131)
siax3 = sifig.add_subplot(132)
siax4 = sifig.add_subplot(133)

# 3. plot the difference between the scenarios with/without UV
ax2.plot(uvdiff.ES_UV, list(reversed(range(13))), "-", color="darkviolet")
ax2.plot(uvdiff.ES_NOUV, list(reversed(range(13))), "k:")
siax2.plot(uvdiff.ES_UV_CYTO, list(reversed(range(13))), "-", color="darkviolet")
siax2.plot(uvdiff.ES_NOUV_CYTO, list(reversed(range(13))), "k:")
for ax in [ax2, siax2]:
    ax.set_ylim(0, 12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("$R_{net}$ (W m$^{-2}$)", fontsize=12)
    ax.set_ylabel("Canopy Layer from Bottom", fontsize=12)
    ax.text(100, 4, "$\Delta$R={:.1f}".format(sum(uvdiff.ES_UV - uvdiff.ES_NOUV)) + " W m$^{-2}$", color="darkviolet", ha="center", va="center")
    ax.set_title("$\\bf{(a)}$", loc="left", fontsize=12)

ax3.plot(uvdiff.TS_UV - 273.15, list(reversed(range(13))), "-", color="darkviolet")
ax3.plot(uvdiff.TS_NOUV - 273.15, list(reversed(range(13))), "k:")
siax3.plot(uvdiff.TS_UV_CYTO - 273.15, list(reversed(range(13))), "-", color="darkviolet")
siax3.plot(uvdiff.TS_NOUV_CYTO - 273.15, list(reversed(range(13))), "k:")
for ax in [ax3, siax3]:
    ax.set_ylim(0, 12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Temperature (℃)", fontsize=12)
    ax.set_title("$\\bf{(b)}$", loc="left", fontsize=12)

ax4.plot(uvdiff.AS_UV, list(reversed(range(13))), "-", color="darkviolet", label="UV")
ax4.plot(uvdiff.AS_NOUV, list(reversed(range(13))), "k:", label="No UV")
siax4.plot(uvdiff.AS_UV_CYTO, list(reversed(range(13))), "-", color="darkviolet", label="UV")
siax4.plot(uvdiff.AS_NOUV_CYTO, list(reversed(range(13))), "k:", label="No UV")
for ax in [ax4, siax4]:
    ax.set_ylim(0, 12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("$A_{g}$ (µmol m$^{-2}$ s$^{-1})$", fontsize=12)
    ax.set_title("$\\bf{(c)}$", loc="left", fontsize=12)
    ax.legend(loc="lower right")

# 4. save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/3_uv_feature.pdf", bbox_inches="tight")
sifig.set_tight_layout(True)
sifig.savefig("../figure/s2_uv_feature.pdf", bbox_inches="tight")
