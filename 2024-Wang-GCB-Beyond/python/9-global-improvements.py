#
# this script is meant to plot the global improvements
#
import matplotlib.pyplot as PLT
import netCDF4 as NC
import numpy as NP
import scipy.stats as STATS


# 0. define the RMSE
def rmse(obs, mod):
    return NP.sqrt(NP.mean((obs - mod) ** 2))

# 1. read the data
dset = NC.Dataset("../output/9_global.nc")
lat = dset.variables["lat"][:]
lon = dset.variables["lon"][:]
gpp_mpi_1y = dset.variables["GPP_MPI_1Y"][:]
gpp_epar_1y = dset.variables["GPP_EPAR_1Y"][:]
gpp_rpar_1y = dset.variables["GPP_RPAR_1Y"][:]
dset.close()

# 2. plot the hexbin comparison between 1Y datasets (epar,rpar ~ mpi)
fig = PLT.figure("9-improve", dpi=300, figsize=(7.5,3.5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
mask = (gpp_epar_1y.data.flatten() > 0.2) & (gpp_rpar_1y.data.flatten() > 0.2) & (gpp_mpi_1y.data.flatten() > 0)
cm1 = ax1.hexbin(gpp_mpi_1y.data.flatten()[mask], gpp_rpar_1y.data.flatten()[mask], gridsize=50, bins="log", cmap="Greys", mincnt=1)
cm2 = ax2.hexbin(gpp_mpi_1y.data.flatten()[mask], gpp_epar_1y.data.flatten()[mask], gridsize=50, bins="log", cmap="Greys", mincnt=1)
fig.colorbar(cm1, ax=ax1)
fig.colorbar(cm2, ax=ax2)

# 3. plot the linear regression
lr1 = STATS.linregress(gpp_mpi_1y.data.flatten()[mask], gpp_rpar_1y.data.flatten()[mask])
lr2 = STATS.linregress(gpp_mpi_1y.data.flatten()[mask], gpp_epar_1y.data.flatten()[mask])
print(lr1)
print(lr2)
x = NP.linspace(0, 10, 100)
ax1.plot(x, lr1.slope*x + lr1.intercept, "r-")
ax2.plot(x, lr2.slope*x + lr2.intercept, "r-")
ax1.text(0, 9, "y = %.3f x + %.3f\nR$^2$ = %.3f\nRMSE = %.3f" % (lr1.slope, lr1.intercept, lr1.rvalue ** 2, rmse(gpp_mpi_1y.data.flatten()[mask], gpp_rpar_1y.data.flatten()[mask])))
ax2.text(0, 9, "y = %.3f x + %.3f\nR$^2$ = %.3f\nRMSE = %.3f" % (lr2.slope, lr2.intercept, lr2.rvalue ** 2, rmse(gpp_mpi_1y.data.flatten()[mask], gpp_epar_1y.data.flatten()[mask])))
ax1.set_ylim(-0.2, 12)
ax2.set_ylim(-0.2, 12)
ax1.set_xlim(-0.2, 10)
ax2.set_xlim(-0.2, 10)
ax1.plot([0, 10], [0, 10], "k:")
ax2.plot([0, 10], [0, 10], "k:")
ax1.set_title("$\\bf{(a)}$ without FR", fontsize=12, loc="left")
ax2.set_title("$\\bf{(b)}$ with FR", fontsize=12, loc="left")
ax1.set_xlabel("MPI RS GPP (gC m$^{-2}$ yr$^{-1}$)", fontsize=12)
ax2.set_xlabel("MPI RS GPP (gC m$^{-2}$ yr$^{-1}$)", fontsize=12)
ax1.set_ylabel("CliMA GPP (gC m$^{-2}$ yr$^{-1}$)", fontsize=12)
"""
LinregressResult(slope=0.847364393341889, intercept=0.46914671666322727, rvalue=0.7920646096097822, pvalue=0.0, stderr=0.005768877730133157, intercept_stderr=0.018219684829213843)
LinregressResult(slope=0.9012928396846578, intercept=0.42617337754956397, rvalue=0.7934336715638266, pvalue=0.0, stderr=0.00610756943991048, intercept_stderr=0.019289365362427008)
"""

# save the figure
fig.set_tight_layout(True)
fig.savefig("../figure/9_global.pdf", bbox_inches="tight")
