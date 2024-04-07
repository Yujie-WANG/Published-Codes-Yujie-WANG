#
# this script is meant for plotting the improvements from the empirical corrections
#
import matplotlib.pyplot as PLT
import pandas as PD


# 1. read the data
df = PD.read_csv("../output/7_improvement.csv")

# 2. create figure 7 with 1 row and 2 columns
fig = PLT.figure(7, dpi=300, figsize=(7.5,3.7))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# 3. plot the GPP improvements for plant 1
ax1.plot(12, 0, "k-", alpha=0.5, label="HS")
ax1.plot(12, 0, "k--", alpha=0.5, label="BB")
ax1.plot(12, 0, "k:", label="BB-corr")
ax1.plot(df.Hour, df.GPP_DEF_1, "c-", alpha=0.5)
ax1.plot(df.Hour, df.GPP_CLM_1, "c--", alpha=0.5)
ax1.plot(df.Hour, df.GPP_EMP_1, "c:")
ax1.set_xlabel("Hour of day", fontsize=12)
ax1.set_ylabel("GPP (μmol m$^{-2}$ s$^{-1}$)", fontsize=12, color="c")
ax1.set_xlim(4, 20)
ax1.set_ylim(0, 42)
ax1.legend()
ax1.set_title("(a)", fontsize=12, loc="left")

# 4. plot the SIF improvements for plant 1
tx1 = ax1.twinx()
tx1.plot(df.Hour, df.SIF_DEF_1, "r-", alpha=0.5)
tx1.plot(df.Hour, df.SIF_CLM_1, "r--", alpha=0.5)
tx1.plot(df.Hour, df.SIF_EMP_1, "r:")
tx1.set_ylim(0, 3)

# 5. plot the improvements for plant 2
ax2.plot(df.Hour, df.GPP_DEF_2, "c-")
ax2.plot(df.Hour, df.GPP_CLM_2, "c--", alpha=0.5)
ax2.plot(df.Hour, df.GPP_EMP_2, "c:")
ax2.set_xlabel("Hour of day", fontsize=12)
ax2.set_xlim(4, 20)
ax2.set_ylim(0, 27)
ax2.set_title("(b)", fontsize=12, loc="left")

# 6. plot the SIF improvements for plant 2
tx2 = ax2.twinx()
tx2.plot(df.Hour, df.SIF_DEF_2, "r-", alpha=0.5)
tx2.plot(df.Hour, df.SIF_CLM_2, "r--", alpha=0.5)
tx2.plot(df.Hour, df.SIF_EMP_2, "r:")
tx2.set_ylim(0, 3)
tx2.set_ylabel("SIF (mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$)", fontsize=12, color="r")

# 7. save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/7_improvement.pdf", bbox_inches="tight")
