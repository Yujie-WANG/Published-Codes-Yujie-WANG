#
# this file is meant to plot the sensitivity of f_apar and f_ppar to the chl, car, lma, lai, and ci
#
import matplotlib.pyplot as PLT
import pandas as PD


# 1. read the data
sensitivity = PD.read_csv("../output/2_sensitivity.csv")

# 2. create figure 2 with 3 rows and 5 columns
fig = PLT.figure(2, dpi=300, figsize=(7.5, 7.5))

# 3. the first row is for f_apar
ax11 = fig.add_subplot(3,5,1)
ax11.plot(sensitivity.FAPAR, list(reversed(range(1,13))), "k:")
ax11.plot(sensitivity.FAPAR_CHL, list(reversed(range(1,13))), "k-", label="Chl ↑")
ax11.set_title("Chl ↑", fontsize=12, loc="center")

ax12 = fig.add_subplot(3,5,2)
ax12.plot(sensitivity.FAPAR, list(reversed(range(1,13))), "k:")
ax12.plot(sensitivity.FAPAR_CAR, list(reversed(range(1,13))), "k-", label="Car ↑")
ax12.set_title("Car ↑", fontsize=12, loc="center")

ax13 = fig.add_subplot(3,5,3)
ax13.plot(sensitivity.FAPAR, list(reversed(range(1,13))), "k:")
ax13.plot(sensitivity.FAPAR_LMA, list(reversed(range(1,13))), "k-", label="LMA ↑")
ax13.set_title("LMA ↑", fontsize=12, loc="center")
ax13.set_xlabel("$f_\mathrm{APAR}$ (-)", fontsize=12)

ax14 = fig.add_subplot(3,5,4)
ax14.plot(sensitivity.FAPAR, list(reversed(range(1,13))), "k:")
ax14.plot(sensitivity.FAPAR_LAI, list(reversed(range(1,13))), "k-", label="LAI ↑")
ax14.set_title("LAI ↑", fontsize=12, loc="center")

ax15 = fig.add_subplot(3,5,5)
ax15.plot(sensitivity.FAPAR, list(reversed(range(1,13))), "k:")
ax15.plot(sensitivity.FAPAR_CI, list(reversed(range(1,13))), "k-", label="CI ↓")
ax15.set_title("CI ↓", fontsize=12, loc="center")

# 4. the second row is for f_ppar
ax21 = fig.add_subplot(3,5,6)
ax21.plot(sensitivity.FPPAR, list(reversed(range(1,13))), "k:")
ax21.plot(sensitivity.FPPAR_CHL, list(reversed(range(1,13))), "k-", label="Chl ↑")

ax22 = fig.add_subplot(3,5,7)
ax22.plot(sensitivity.FPPAR, list(reversed(range(1,13))), "k:")
ax22.plot(sensitivity.FPPAR_CAR, list(reversed(range(1,13))), "k-", label="Car ↑")

ax23 = fig.add_subplot(3,5,8)
ax23.plot(sensitivity.FPPAR, list(reversed(range(1,13))), "k:")
ax23.plot(sensitivity.FPPAR_LMA, list(reversed(range(1,13))), "k-", label="LMA ↑")
ax23.set_xlabel("$f_\mathrm{PPAR}$ (-)", fontsize=12)

ax24 = fig.add_subplot(3,5,9)
ax24.plot(sensitivity.FPPAR, list(reversed(range(1,13))), "k:")
ax24.plot(sensitivity.FPPAR_LAI, list(reversed(range(1,13))), "k-", label="LAI ↑")

ax25 = fig.add_subplot(3,5,10)
ax25.plot(sensitivity.FPPAR, list(reversed(range(1,13))), "k:")
ax25.plot(sensitivity.FPPAR_CI, list(reversed(range(1,13))), "k-", label="CI ↓")

# 5. the third row is for f_apar * f_ppar
ax31 = fig.add_subplot(3,5,11)
ax31.plot(sensitivity.FAPAR * sensitivity.FPPAR, list(reversed(range(1,13))), "k:")
ax31.plot(sensitivity.FAPAR_CHL * sensitivity.FPPAR_CHL, list(reversed(range(1,13))), "k-", label="Chl ↑")

ax32 = fig.add_subplot(3,5,12)
ax32.plot(sensitivity.FAPAR * sensitivity.FPPAR, list(reversed(range(1,13))), "k:")
ax32.plot(sensitivity.FAPAR_CAR * sensitivity.FPPAR_CAR, list(reversed(range(1,13))), "k-", label="Car ↑")

ax33 = fig.add_subplot(3,5,13)
ax33.plot(sensitivity.FAPAR * sensitivity.FPPAR, list(reversed(range(1,13))), "k:")
ax33.plot(sensitivity.FAPAR_LMA * sensitivity.FPPAR_LMA, list(reversed(range(1,13))), "k-", label="LMA ↑")
ax33.set_xlabel("$f_\mathrm{APAR} \cdot f_\mathrm{PPAR}$ (-)", fontsize=12)

ax34 = fig.add_subplot(3,5,14)
ax34.plot(sensitivity.FAPAR * sensitivity.FPPAR, list(reversed(range(1,13))), "k:")
ax34.plot(sensitivity.FAPAR_LAI * sensitivity.FPPAR_LAI, list(reversed(range(1,13))), "k-", label="LAI ↑")

ax35 = fig.add_subplot(3,5,15)
ax35.plot(sensitivity.FAPAR * sensitivity.FPPAR, list(reversed(range(1,13))), "k:")
ax35.plot(sensitivity.FAPAR_CI * sensitivity.FPPAR_CI, list(reversed(range(1,13))), "k-", label="CI ↓")

# 6. set the limits and labels
for ax in [ax11, ax12, ax13, ax14, ax15, ax21, ax22, ax23, ax24, ax25, ax31, ax32, ax33, ax34, ax35]:
    ax.set_ylim(0.8, 12.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

for ax in [ax12, ax13, ax14, ax15, ax22, ax23, ax24, ax25, ax32, ax33, ax34, ax35]:
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.spines["left"].set_visible(False)

for ax in [ax11, ax12, ax13, ax14, ax15]:
    ax.set_xlim(0.55,0.85)
    ax.set_xticks([0.6, 0.7, 0.8])

for ax in [ax21, ax22, ax23, ax24, ax25]:
    ax.set_xlim(0.71, 0.86)

for ax in [ax31, ax32, ax33, ax34, ax35]:
    ax.set_xlim(0.44, 0.73)
    ax.set_xticks([0.5, 0.6, 0.7])

for ax in [ax11, ax21, ax31]:
    ax.set_yticks(range(1,13))
    ax.set_yticklabels(["", "2", "", "4", "", "6", "", "8", "", "10", "", "12"])
    ax.set_ylabel("Layer from bottom", fontsize=12)

# 7. save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/2_sensitivity.pdf", bbox_inches="tight")
