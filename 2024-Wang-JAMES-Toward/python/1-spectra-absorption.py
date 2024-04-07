#
# this script will plot the absorption spectra of the pigments and leaf absorptance with different chlorophyll content
#
import matplotlib.pyplot as PLT
import pandas as PD


# 1. read the data
spectra = PD.read_csv("../output/1_spectra.csv")
chl_abs = PD.read_csv("../output/1_absorption.csv")
mask_par = (spectra.WL >= 400) & (spectra.WL <= 700)

# 2. create figure 1 to plot on
fig = PLT.figure(1, dpi=300, figsize=(7.5, 3.5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# 3. plot the absorption spectra of the pigments
ax1.plot(spectra.WL[mask_par], spectra.KCHL[mask_par], "-", label="Chl", color="forestgreen")
ax1.plot(spectra.WL[mask_par], spectra.KCARV[mask_par], "-", label="Car", color="orange")
ax1.plot(spectra.WL[mask_par], spectra.KCARZ[mask_par], "-", color="orange")
ax1.plot(spectra.WL[mask_par], spectra.KH2O[mask_par], "-", label="H₂O", color="c")
ax1.plot(spectra.WL[mask_par], spectra.KLMA[mask_par] / 1000, "-", label="LMA", color="brown")
ax1.legend()
ax1.set_xlim(400, 700)
ax1.set_xlabel("Wavelength (nm)", fontsize=12)
ax1.set_ylabel("Absorption Coefficient (-)", fontsize=12)
ax1.set_title("(a)", fontsize=12, loc="left")

# 4. plot the leaf absorptance with different chlorophyll content
ax2.plot(chl_abs.CHL, chl_abs.F_APAR_NATURE * chl_abs.F_PPAR_NATURE, "-", label="Nature", color="k")
ax2.plot(chl_abs.CHL, chl_abs.F_APAR_BLUE * chl_abs.F_PPAR_BLUE, "-", label="Blue", color="b", alpha=0.7)
ax2.plot(chl_abs.CHL, chl_abs.F_APAR_GREEN * chl_abs.F_PPAR_GREEN, "-", label="Green", color="g", alpha=0.7)
ax2.plot(chl_abs.CHL, chl_abs.F_APAR_RED * chl_abs.F_PPAR_RED, "-", label="Red", color="r", alpha=0.7)
ax2.legend()
ax2.set_xlim(0, 50)
ax2.set_ylim(0, 1)
ax2.set_xlabel("Chlorophyll Content (µg cm$^{-2}$)", fontsize=12)
ax2.set_ylabel("$f_\mathrm{APAR} \cdot f_\mathrm{PPAR}$ (-)", fontsize=12)
ax2.set_title("(b)", fontsize=12, loc="left")

# save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/1_spectra.pdf", bbox_inches="tight")
