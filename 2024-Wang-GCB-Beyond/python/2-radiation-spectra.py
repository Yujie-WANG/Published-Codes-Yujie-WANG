#
# this script is used to plot the radiation spectra for different layers of a canopy
#
import matplotlib.pyplot as PLT
import pandas as PD


# 1. read the data
data = PD.read_csv("../output/2_radiation_spectra.csv")

# 2. create figure 3
fig = PLT.figure(3, figsize=(6, 6))
ax = fig.add_subplot(111)

# 3. plot the data
for i in range(1, 13):
    ax.plot(data["WL"], data["SPECTRUM_{}".format(i)], "k-", alpha=(16-i)/16, linewidth=1)
ax.fill_between([700,750], 0, 1800, color="red", alpha=0.2)
ax.fill_between([510,560], 0, 1800, color="green", alpha=0.2)
ax.set_xlim(300, 800)
ax.set_ylim(0, 1800)
ax.set_xlabel("Wavelength (nm)", fontsize=12)
ax.set_ylabel("Radiation (mW m$^{-2}$ nm$^{-1}$)", fontsize=12)

# 4. save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/2_radiation_spectra.pdf", bbox_inches="tight")
