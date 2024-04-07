#
# this script is meant for preparing the data for a SI figure to show the seasonality of ΔGPP and ΔSIF
#
#
# this file is meant for compute the global scale changes in GPP and SIF to plot in Python
#
using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: read_LUT
using NetcdfIO: append_nc!, read_nc, save_nc!

using Emerald.EmeraldMath.Stats: nanmean


# 1. define where the files are located
clima_land_sims = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
sim_tag = "a11_gm3_wd1";
hs_kdtd = "$(clima_land_sims)/$(sim_tag)_2019_1X_1M.hs.c3.epar.kdtd.nc";
bb_kdst = "$(clima_land_sims)/$(sim_tag)_2019_1X_1M.bb.c3.epar.kdst.nc";
lndmask = regrid(read_LUT(query_collection("LM_4X_1Y_V1"))[1], 1);

# 2. read the data and mask the NaNs to 0 if land mask > 0.5
gpp_hs_kdtd = read_nc(hs_kdtd, "mGPP") .* (3600 * 24 * 1e-6 * 12);
gpp_bb_kdst = read_nc(bb_kdst, "mGPP") .* (3600 * 24 * 1e-6 * 12);
sif_hs_kdtd = read_nc(hs_kdtd, "mSIF740");
sif_bb_kdst = read_nc(bb_kdst, "mSIF740");
δgpp = gpp_bb_kdst .- gpp_hs_kdtd;
δsif = sif_bb_kdst .- sif_hs_kdtd;
nan_mask = isnan.(gpp_hs_kdtd[:,:,1]) .&& (lndmask .> 0.1);

# 3. compute the seasonal mean difference
csons = [[3,4,5], [6,7,8], [9,10,11], [12,1,2]];
gpp_maps = zeros(360,180,4);
sif_maps = zeros(360,180,4);
for i in eachindex(csons)
    gpp_map = ones(360,180);
    sif_map = ones(360,180);
    for ilon in 1:360, ilat in 1:180
        gpp_map[ilon,ilat] = nanmean(δgpp[ilon,ilat,csons[i]]);
        sif_map[ilon,ilat] = nanmean(δsif[ilon,ilat,csons[i]]);
    end;
    gpp_map[nan_mask] .= 0;
    sif_map[nan_mask] .= 0;
    gpp_maps[:,:,i] .= gpp_map;
    sif_maps[:,:,i] .= sif_map;
end;

# 4. save the data
save_nc!("../output/s1_seasonality.nc", "ΔGPP", gpp_maps, Dict{String,String}("about" => "Broadband mode with constant KD GPP - Hyperspectral mode with KDTD GPP"));
append_nc!("../output/s1_seasonality.nc", "ΔSIF", sif_maps, Dict{String,String}("about" => "Broadband mode with constant KD SIF - Hyperspectral mode with KDTD SIF"), ["lon", "lat", "ind"]);
