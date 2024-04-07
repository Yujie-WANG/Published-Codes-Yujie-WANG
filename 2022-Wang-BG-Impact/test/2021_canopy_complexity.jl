# information of time
using PkgUtility

@info tinfo("Running the script...");

# use this to disable the matplotlib warning on server
ENV["MPLBACKEND"] = "Agg";

using ResearchProjects

FT   = Float64;
proj = CanopyComplexity2021{FT}();
plot_CC2021!(proj; saving=true, use_latex=true);

# information of time
@info tinfo("Finished running the script!");
