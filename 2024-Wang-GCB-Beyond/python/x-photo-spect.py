#
# this file is meant to plot the photo spect related simulations
#
import matplotlib.pyplot as PLT
import netCDF4 as NC
import numpy as NP


# 1. read the data, shape of the data is
#     0: wavelength
#     1: time
dset = NC.Dataset("../output/6_ps_spectra.nc", "r")
nthfile = dset.variables["file"][:]
time = dset.variables["time"][:]
wl = dset.variables["wl"][:]
nat_spectra = dset.variables["spectra"][:]
dset.close()
dwl = wl[1:len(wl)] - wl[0:len(wl)-1]

"""
# 2. plot the spectra per day
for i in range(max(nthfile)):
    print("Plotting the spectra for %dth day..." % i)
    subset = nat_spectra[:,nthfile==i+1]
    fig = PLT.figure(i, dpi=300, figsize=(6,4))
    ax = fig.add_subplot(111)
    for j in range(subset.shape[1]):
        ax.plot(wl, subset[:,j], "k-", alpha=0.3, linewidth=0.5)
    fig.savefig("../output/test_%03d.png" % i, bbox_inches="tight")

# 2. plot the data on a figure
#     ax1: mean spectra
#     ax2: total radiation
fig = PLT.figure("X", dpi=300, figsize=(12,4))
ax1 = fig.add_subplot(121)
mean_spectra = []
for i in range(nat_spectra.shape[0]):
    mean_spectra.append(NP.mean(nat_spectra[:,i]))

print("Plotting the mean spectra")
ax1.plot(wl, mean_spectra, "k-")

ax2 = fig.add_subplot(122)
total_rad = []
for i in range(nat_spectra.shape[1]):
    mean_rd = (nat_spectra[0:len(wl)-1,i] + nat_spectra[1:len(wl),i]) / 2
    total_rad.append(sum(mean_rd * dwl))

print("Plotting the total radiation")
ax2.hexbin(time % 1, total_rad, gridsize=100, cmap="Greys")
# ax2.plot(time % 1, total_rad, "k.", alpha=0.2)

fig.set_tight_layout(True)
fig.savefig("../output/test.png", bbox_inches="tight")
"""
