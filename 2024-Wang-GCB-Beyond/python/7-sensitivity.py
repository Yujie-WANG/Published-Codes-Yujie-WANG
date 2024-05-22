#
# this script is meant to plot the sensitivity analysis of the model
#
import matplotlib.pyplot as PLT
import pandas as PD


# 1. read the data
df = PD.read_csv("../output/7_sensitivity.csv")

# 2. create figure 6 with 1 row and 4 columns
fig = PLT.figure(6, dpi=300, figsize=(12,3.5))
ax1 = fig.add_subplot(141)
ax2 = fig.add_subplot(142)
ax3 = fig.add_subplot(143)
ax4 = fig.add_subplot(144)

# 3. plot the data
ax1.plot(df.CHL, df.GPP_750_CHL - df.GPP_700_CHL, "k-")
ax2.plot(df.LAI, df.GPP_750_LAI - df.GPP_700_LAI, "k-")
ax3.plot(df.CI, df.GPP_750_CI - df.GPP_700_CI, "k-")
ax4.plot(df.SZA, df.GPP_750_SZA - df.GPP_700_SZA, "k-")
ax4.plot(df.SZA, df.GPP_750_SZA_2 - df.GPP_700_SZA_2, "k:")
ax1.set_ylabel("$\Delta$GPP (μmol m$^{-2}$ s$^{-1}$)", fontsize=12)
ax1.set_xlabel("Chlorophyll content (μg cm$^{-2}$)", fontsize=12)
ax2.set_xlabel("Leaf area index (-)", fontsize=12)
ax3.set_xlabel("Clumping index (-)", fontsize=12)
ax4.set_xlabel("Solar zenith angle ($^{\circ}$)", fontsize=12)
ax1.set_title("$\\bf{(a)}$", fontsize=12, loc="left")
ax2.set_title("$\\bf{(b)}$", fontsize=12, loc="left")
ax3.set_title("$\\bf{(c)}$", fontsize=12, loc="left")
ax4.set_title("$\\bf{(d)}$", fontsize=12, loc="left")

# 4. save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/7_sensitivity.pdf", bbox_inches="tight")
