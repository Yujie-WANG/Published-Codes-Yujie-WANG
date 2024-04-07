#
# this file is meant for compute the global scale changes in GPP and SIF to plot in Python
#
using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: read_LUT
using NetcdfIO: append_nc!, read_nc, save_nc!


# 1. define where the files are located
clima_land_sims = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
sim_tag = "a11_gm3_wd1";
hs_kdtd = "$(clima_land_sims)/$(sim_tag)_2019_1X_1Y.hs.c3.epar.kdtd.nc";
hs_kdst = "$(clima_land_sims)/$(sim_tag)_2019_1X_1Y.hs.c3.epar.kdst.nc";
bb_kdtd = "$(clima_land_sims)/$(sim_tag)_2019_1X_1Y.bb.c3.epar.kdtd.nc";
bb_kdst = "$(clima_land_sims)/$(sim_tag)_2019_1X_1Y.bb.c3.epar.kdst.nc";
lndmask = regrid(read_LUT(query_collection("LM_4X_1Y_V1"))[1], 1);

# 2. read the data and mask the NaNs to 0 if land mask > 0.5
gpp_hs_kdtd = read_nc(hs_kdtd, "mGPP") .* (3600 * 24 * 1e-6 * 12);
gpp_hs_kdst = read_nc(hs_kdst, "mGPP") .* (3600 * 24 * 1e-6 * 12);
gpp_bb_kdtd = read_nc(bb_kdtd, "mGPP") .* (3600 * 24 * 1e-6 * 12);
gpp_bb_kdst = read_nc(bb_kdst, "mGPP") .* (3600 * 24 * 1e-6 * 12);
sif_hs_kdtd = read_nc(hs_kdtd, "mSIF740");
sif_hs_kdst = read_nc(hs_kdst, "mSIF740");
sif_bb_kdtd = read_nc(bb_kdtd, "mSIF740");
sif_bb_kdst = read_nc(bb_kdst, "mSIF740");
mask = (lndmask .> 0.5) .&& isnan.(gpp_hs_kdtd);
gpp_hs_kdtd[mask] .= 0;
gpp_hs_kdst[mask] .= 0;
gpp_bb_kdtd[mask] .= 0;
gpp_bb_kdst[mask] .= 0;
sif_hs_kdtd[mask] .= 0;
sif_hs_kdst[mask] .= 0;
sif_bb_kdtd[mask] .= 0;
sif_bb_kdst[mask] .= 0;

# 3. save the data
save_nc!("../output/4_gpp_sif.nc", "GPP_HS_KDTD", gpp_hs_kdtd, Dict{String,String}("about" => "GPP - Hyperspectral mode with KD TD"));
append_nc!("../output/4_gpp_sif.nc", "GPP_HS_KDST", gpp_hs_kdst, Dict{String,String}("about" => "GPP - Hyperspectral mode with constant KD"), ["lon", "lat"]);
append_nc!("../output/4_gpp_sif.nc", "GPP_BB_KDTD", gpp_bb_kdtd, Dict{String,String}("about" => "GPP - Broadband mode with KD TD"), ["lon", "lat"]);
append_nc!("../output/4_gpp_sif.nc", "GPP_BB_KDST", gpp_bb_kdst, Dict{String,String}("about" => "GPP - Broadband mode with constant KD"), ["lon", "lat"]);
append_nc!("../output/4_gpp_sif.nc", "SIF_HS_KDTD", sif_hs_kdtd, Dict{String,String}("about" => "SIF - Hyperspectral mode with KD TD"), ["lon", "lat"]);
append_nc!("../output/4_gpp_sif.nc", "SIF_HS_KDST", sif_hs_kdst, Dict{String,String}("about" => "SIF - Hyperspectral mode with constant KD"), ["lon", "lat"]);
append_nc!("../output/4_gpp_sif.nc", "SIF_BB_KDTD", sif_bb_kdtd, Dict{String,String}("about" => "SIF - Broadband mode with KD TD"), ["lon", "lat"]);
append_nc!("../output/4_gpp_sif.nc", "SIF_BB_KDST", sif_bb_kdst, Dict{String,String}("about" => "SIF - Broadband mode with constant KD"), ["lon", "lat"]);
