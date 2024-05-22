#
# This script is meant to plot the model output to illustrate that the new model is working as expected
#
import pandas as PD
import matplotlib.pyplot as PLT

# 1. read the data
nr = PD.read_csv("../output/1_spectra.csv")
fitting = PD.read_csv("../output/1_new_spectra.csv")

# 2. create a figure to plot on
fig = PLT.figure(1, dpi=300, figsize=(6,6))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

# 3. plot the natural radiation on ax1
ax1.plot(nr.WL, nr.E_DIR, 'k-', label="direct")
ax1.plot(nr.WL, nr.E_DIF, 'k:', label="diffuse")
ax1.fill_between([300,400], 0, 1450, color="darkviolet", alpha=0.2)
ax1.fill_between([400,700], 0, 1450, color="green", alpha=0.2)
ax1.fill_between([700,750], 0, 1450, color="red", alpha=0.2)
rad_uv = (sum(nr.E_DIR[nr.WL < 400]) + sum(nr.E_DIF[nr.WL < 400])) / 1000
rad_pr = (sum(nr.E_DIR[(nr.WL >= 400) & (nr.WL < 700)]) + sum(nr.E_DIF[(nr.WL >= 400) & (nr.WL < 700)])) / 1000
rad_fr = (sum(nr.E_DIR[(nr.WL >= 700) & (nr.WL < 750)]) + sum(nr.E_DIF[(nr.WL >= 700) & (nr.WL < 750)])) / 1000
pho_uv = 138.33208517312244
pho_pr = 1983.5032262756463
pho_fr = 363.93079663899493
ax1.text(360, 800, "UV\n{:.1f} W m⁻²\n{:.1f} μmol m⁻²".format(rad_uv,pho_uv), fontsize=12, ha="center", va="center", color="darkviolet")
ax1.text(520, 640, "PAR\n{:.1f} W m⁻²\n{:.1f} μmol m⁻²".format(rad_pr,pho_pr), fontsize=12, ha="center", va="center", color="green")
ax1.text(700, 640, "FR\n{:.1f} W m⁻²\n{:.1f} μmol m⁻²".format(rad_fr,pho_fr), fontsize=12, ha="center", va="center", color="red")
ax1.legend(loc="best")
ax1.set_xlim(300, 750)
ax1.set_ylim(0, 1450)
ax1.set_ylabel("Flux (mW m$^{-2}$ nm$^{-1}$)", fontsize=12)
ax1.set_title("$\\bf{(a)}$", fontsize=12, loc="left")

# 4. plot the absorption features on ax2
mask = (400 <= nr.WL) & (nr.WL <= 750)
mask_uv = fitting.WL <= 400
ax2.plot(nr.WL[mask], nr.REF_CHL[mask] * 40, "-", color="forestgreen", label="Chlorophyll")
ax2.plot(nr.WL[mask], nr.REF_CARV[mask] * 40 / 7, "-", color="orange", label="Carotenoid")
ax2.plot(nr.WL[mask], nr.REF_CARZ[mask] * 40 / 7, "-", color="orange")
ax2.plot(nr.WL[mask], nr.REF_LMA[mask] * 0.012, "-", color="brown", label="Dry matter")
ax2.plot(nr.WL[mask], nr.REF_WATER[mask] * 5 * 18 / 10000, "-", color="c", label="Water")
ax2.plot(400, 0, "k:", label="Reconstructed")
ax2.plot(fitting.WL[mask_uv], fitting.REF_CHL[mask_uv] * 40, ":", color="forestgreen")
ax2.plot(fitting.WL[mask_uv], fitting.REF_CARV[mask_uv] * 40 / 7, ":", color="orange")
ax2.plot(fitting.WL[mask_uv], fitting.REF_CARZ[mask_uv] * 40 / 7, ":", color="orange")
ax2.plot(fitting.WL[mask_uv], fitting.REF_LMA[mask_uv] * 0.012, ":", color="brown")
ax2.plot(fitting.WL[mask_uv], fitting.REF_H2O[mask_uv] * 5 * 18 / 10000, ":", color="c")
ax2.fill_between([300,400], -0.05, 3.1, color="darkviolet", alpha=0.2)
ax2.fill_between([700,750], -0.05, 3.1, color="red", alpha=0.2)
ax2.legend(loc="best")
ax2.set_xlim(300, 750)
ax2.set_ylim(-0.05, 3.1)
ax2.set_xlabel("Wavelength (nm)", fontsize=12)
ax2.set_ylabel("Relative Absorption Feature (-)", fontsize=12)
ax2.set_title("$\\bf{(b)}$", fontsize=12, loc="left")

# 5. save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/1_uv_fr_problem.pdf", bbox_inches="tight")
