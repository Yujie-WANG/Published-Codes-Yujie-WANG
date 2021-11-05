# information of time
using PkgUtility

@info tinfo("Running the script...");

# use this to disable the matplotlib warning on server
ENV["MPLBACKEND"] = "Agg";

using ResearchProjects

# define general data, use constant g1 for empirical model
FT    = Float64;
proj  = ForestTower2021{FT}();
rerun = false;
NTH   = 130;

# rerun the fittings, manually run it on the server
if rerun
    # run this on server curry
    fit_tower!(proj, NTH);
    # copy the results to local PC and run these on local PC
    divide_results!(proj);
    sif_series!(proj);
end

# plot the results with fitted Gmax
plot_FT2021!(proj; saving=true);
plot_FT2021_SI!(proj; saving=true);

# information of time
@info tinfo("Finished running the script!");
