#
# this script is meant to plot the comparision of the reflected SW radiation
#
import matplotlib.pyplot as PLT
import numpy as NP
import pandas as PD
import scipy.stats as STATS


# 0. define the RMSE
def rmse(obs, mod):
    return NP.sqrt(NP.mean((obs - mod) ** 2))

# 1. read the data from the file 11_flux_tower.csv
df = PD.read_csv("../output/4_flux_tower.csv")
times = (NP.array(list(range(31*48))) + 0.5) * 0.5 / 24
mask = (df.SW_OUT_UV >= 0) & (df.SW_OUT_NOUV >= 0)

# 2. calculate the half-hourly means
hours = (NP.array(list(range(48))) + 0.5) * 0.5
means_obs = []
means_UV = []
means_NOUV = []
std_obs = []
std_UV = []
std_NOUV = []
for i in range(48):
    dat_obs = df.SW_OUT[i:len(df):48]
    dat_UV = df.SW_OUT_UV[i:len(df):48]
    dat_NOUV = df.SW_OUT_NOUV[i:len(df):48]
    mask_dat = (dat_UV >= 0) & (dat_NOUV >= 0)
    means_obs.append(NP.mean(dat_obs[mask_dat]))
    means_UV.append(NP.mean(dat_UV[mask_dat]))
    means_NOUV.append(NP.mean(dat_NOUV[mask_dat]))
    std_obs.append(NP.std(dat_obs[mask_dat]))
    std_UV.append(NP.std(dat_UV[mask_dat]))
    std_NOUV.append(NP.std(dat_NOUV[mask_dat]))

# 3. plot the data
fig = PLT.figure(11, dpi=300, figsize=(12.5, 6.5))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(223)
ax3 = fig.add_subplot(247)
ax4 = fig.add_subplot(248)

ax1.plot(times[mask], df.SW_OUT[mask], "k-", alpha=0.5, label="obs")
ax1.plot(times[mask], df.SW_OUT_UV[mask], "r-", alpha=0.5, label="UV")
ax1.plot(times[mask], df.SW_OUT_NOUV[mask], "b-", alpha=0.5, label="no UV")
ax1.set_xlim(0, 31)
ax1.set_xlabel("Day of July 2019", fontsize=12)
ax1.set_ylabel("SW out (W m$^{-2}$)", fontsize=12)
ax1.set_title("$\\bf{(a)}$", fontsize=12, loc="left")

ax2.plot(hours, means_obs, "k-", alpha=0.5, label="obs")
ax2.plot(hours, means_UV, "r-", alpha=0.5, label="UV")
ax2.plot(hours, means_NOUV, "b-", alpha=0.5, label="no UV")
ax2.fill_between(hours, NP.array(means_obs) - NP.array(std_obs), NP.array(means_obs) + NP.array(std_obs), color="k", alpha=0.2)
ax2.fill_between(hours, NP.array(means_UV) - NP.array(std_UV), NP.array(means_UV) + NP.array(std_UV), color="r", alpha=0.2)
ax2.fill_between(hours, NP.array(means_NOUV) - NP.array(std_NOUV), NP.array(means_NOUV) + NP.array(std_NOUV), color="b", alpha=0.2)
ax2.legend(loc="upper right")
ax2.set_xlim(0, 24)
ax2.set_xlabel("Hour of day", fontsize=12)
ax2.set_ylabel("SW out (W m$^{-2}$)", fontsize=12)
ax2.set_title("$\\bf{(b)}$", fontsize=12, loc="left")

mask = (df.SW_OUT_UV > 0) & (df.SW_OUT_NOUV > 0)

ax3.plot(df.SW_OUT[mask], df.SW_OUT_UV[mask], "k.", alpha=0.2, mfc="none")
ax3.plot([0,110], [0,110], "k:")
ax3.set_xlim(0, 110)
ax3.set_ylim(0, 110)
ax3.set_xlabel("obs SW out (W m$^{-2}$)", fontsize=12)
ax3.set_ylabel("UV SW out (W m$^{-2}$)", fontsize=12)
ax3.set_title("$\\bf{(c)}$", fontsize=12, loc="left")
lr3 = STATS.linregress(df.SW_OUT[mask], df.SW_OUT_UV[mask])
ax3.plot([0,110], [lr3.intercept, lr3.slope*110 + lr3.intercept], "r-")
ax3.text(5, 85, "y = %.3f x $-$ %.3f\nR$^2$ = %.3f\nRMSE = %.3f" % (lr3.slope, -lr3.intercept, lr3.rvalue ** 2, rmse(df.SW_OUT[mask], df.SW_OUT_UV[mask])))

ax4.plot(df.SW_OUT[mask], df.SW_OUT_NOUV[mask], "k.", alpha=0.2, mfc="none")
ax4.plot([0,110], [0,110], "k:")
ax4.set_xlim(0, 110)
ax4.set_ylim(0, 110)
ax4.set_xlabel("obs SW out (W m$^{-2}$)", fontsize=12)
ax4.set_ylabel("no UV SW out (W m$^{-2}$)", fontsize=12)
ax4.set_title("$\\bf{(d)}$", fontsize=12, loc="left")
lr4 = STATS.linregress(df.SW_OUT[mask], df.SW_OUT_NOUV[mask])
ax4.plot([0,110], [lr4.intercept, lr4.slope*110 + lr4.intercept], "b-")
ax4.text(5, 85, "y = %.3f x $-$ %.3f\nR$^2$ = %.3f\nRMSE = %.3f" % (lr4.slope, -lr4.intercept, lr4.rvalue ** 2, rmse(df.SW_OUT[mask], df.SW_OUT_NOUV[mask])))

# 3. save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/4_albedo.pdf")
