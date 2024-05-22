#
# this script is used to prepare data to plot global simulation benchmarks
#
using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using NetcdfIO: append_nc!, read_nc, save_nc!

using Emerald.EmeraldMath.Regression: linear_regress
using Emerald.EmeraldMath.Stats: nanmean


# 1. define where the files are located
clima_land_sims = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
sim_tag = "a12_gm3_wd1";
hs_epar = "$(clima_land_sims)/$(sim_tag)_2019_1X_1Y.c3.hs.epar.kdtd.oldphi.int.nc";
hs_rpar = "$(clima_land_sims)/$(sim_tag)_2019_1X_1Y.c3.hs.rpar.kdtd.oldphi.int.nc";
lndmask = regrid(read_nc(query_collection("LM_4X_1Y_V1"), "data"), 1);

# 2. read the data
gpp_hs_epar = read_nc(hs_epar, "mGPP") .* (3600 * 24 * 1e-6 * 12);
gpp_hs_rpar = read_nc(hs_rpar, "mGPP") .* (3600 * 24 * 1e-6 * 12);
mask = (lndmask .> 0.1) .&& isnan.(gpp_hs_epar);
gpp_hs_epar[mask] .= 0;
gpp_hs_rpar[mask] .= 0;

# 3. read the data from monthly results
hs_epar_1m = "$(clima_land_sims)/$(sim_tag)_2019_1X_1M.c3.hs.epar.kdtd.oldphi.int.nc";
hs_rpar_1m = "$(clima_land_sims)/$(sim_tag)_2019_1X_1M.c3.hs.rpar.kdtd.oldphi.int.nc";
gpp_hs_epar_1m = read_nc(hs_epar_1m, "mGPP") .* (3600 * 24 * 1e-6 * 12);
gpp_hs_rpar_1m = read_nc(hs_rpar_1m, "mGPP") .* (3600 * 24 * 1e-6 * 12);
for imon in 1:12
    mask_1m = (lndmask .> 0.1) .&& isnan.(gpp_hs_epar_1m[:,:,imon]);
    gpp_hs_epar_1m[mask,imon] .= 0;
    gpp_hs_rpar_1m[mask,imon] .= 0;
end;

# 3. read the global GPP map
gpp_mpi_1m = regrid(read_nc(query_collection("GPP_MPI_RS_2X_1M_2019_V1"), "data"), 1);
gpp_mpi_1y = ones(360,180) .* NaN;
for ilon in 1:360, ilat in 1:180
    gpp_mpi_1y[ilon,ilat] = nanmean(gpp_mpi_1m[ilon,ilat,:]);
end;

# 7. save the data
outfile = "../output/9_global.nc";
if !isfile(outfile)
    save_nc!(outfile, "GPP_MPI_1M", gpp_mpi_1m, Dict{String,String}("about" => "1M GPP from MODIS 2X 1M 2019"));
    append_nc!(outfile, "GPP_EPAR_1M", gpp_hs_epar_1m, Dict{String,String}("about" => "1M GPP from CLIMA-LAND 1X 1M 2019"), ["lon", "lat", "ind"]);
    append_nc!(outfile, "GPP_RPAR_1M", gpp_hs_rpar_1m, Dict{String,String}("about" => "1M GPP from CLIMA-LAND 1X 1M 2019"), ["lon", "lat", "ind"]);
    append_nc!(outfile, "GPP_MPI_1Y", gpp_mpi_1y, Dict{String,String}("about" => "1Y GPP from MODIS 2X 1M 2019"), ["lon", "lat"]);
    append_nc!(outfile, "GPP_EPAR_1Y", gpp_hs_epar, Dict{String,String}("about" => "1Y GPP from CLIMA-LAND 1X 1Y 2019"), ["lon", "lat"]);
    append_nc!(outfile, "GPP_RPAR_1Y", gpp_hs_rpar, Dict{String,String}("about" => "1Y GPP from CLIMA-LAND 1X 1Y 2019"), ["lon", "lat"]);
end;
