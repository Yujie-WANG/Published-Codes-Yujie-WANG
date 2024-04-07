#
# this script is meant for plotting the empirical relationship between modeled and fitted f_ppar
#
import matplotlib.pyplot as PLT
import pandas as PD


# 1. read the data
df = PD.read_csv("../output/6_empirical.csv")

# 2. create figure 6 with 1 row and 2 columns
fig = PLT.figure(6, dpi=300, figsize=(7.5,3.3))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# 3. plot the empirical relationship
cm1 = ax1.hexbin(df.F_PPAR_SL, df.F_PPAR_SL_PRED, gridsize=50, cmap="tab20c_r", mincnt=1)
ax1.plot([0.52, 0.84], [0.52, 0.84], "r:")
ax1.set_xlabel("Modeled $f_\mathrm{PPAR}$", fontsize=12)
ax1.set_ylabel("Fitted $f_\mathrm{PPAR}$", fontsize=12)
ax1.set_xlim(0.52, 0.84)
ax1.set_ylim(0.52, 0.84)
ax1.set_title("(a) sunlit leaves", fontsize=12, loc="left")
fig.colorbar(cm1, ax=ax1, label="count")

cm2 = ax2.hexbin(df.F_PPAR_SH, df.F_PPAR_SH_PRED, gridsize=50, cmap="tab20c_r", mincnt=1)
ax2.plot([0.52, 0.84], [0.52, 0.84], "r:")
ax2.set_xlabel("Modeled $f_\mathrm{PPAR}$", fontsize=12)
ax2.set_xlim(0.52, 0.84)
ax2.set_ylim(0.52, 0.84)
ax2.set_title("(b) shaded leaves", fontsize=12, loc="left")
fig.colorbar(cm2, ax=ax2, label="count")

# 4. save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/6_empirical.pdf", bbox_inches="tight")
