#
# this script is used to prepare data to plot global simulation results
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
dgpp = gpp_hs_epar .- gpp_hs_rpar;
rgpp = (gpp_hs_epar ./ gpp_hs_rpar .- 1) * 100;
mask_r = (lndmask .> 0.1) .&& isnan.(rgpp);
rgpp[mask_r] .= 0;

# 3. read the global LAI map
lai_8d = regrid(read_nc(query_collection("LAI_MODIS_2X_8D_2019_V1"), "data"), 1);
lai_1y = ones(360,180) .* NaN;
for ilon in 1:360, ilat in 1:180
    lai_1y[ilon,ilat] = nanmean(lai_8d[ilon,ilat,:]);
end;

# 4. read the global CI map
ci_1m = regrid(read_nc(query_collection("CI_2X_1M_V3"), "data"), 1);
ci_1y = ones(360,180) .* NaN;
for ilon in 1:360, ilat in 1:180
    ci_1y[ilon,ilat] = nanmean(ci_1m[ilon,ilat,:]);
end;

# 5. read the global CHL map
chl_7d = regrid(read_nc(query_collection("CHL_2X_7D_V1"), "data"), 1);
chl_1y = ones(360,180) .* NaN;
for ilon in 1:360, ilat in 1:180
    chl_1y[ilon,ilat] = nanmean(chl_7d[ilon,ilat,:]);
end;

# 6. read the data from monthly results
hs_epar_1m = "$(clima_land_sims)/$(sim_tag)_2019_1X_1M.c3.hs.epar.kdtd.oldphi.int.nc";
hs_rpar_1m = "$(clima_land_sims)/$(sim_tag)_2019_1X_1M.c3.hs.rpar.kdtd.oldphi.int.nc";
gpp_hs_epar_1m = read_nc(hs_epar_1m, "mGPP") .* (3600 * 24 * 1e-6 * 12);
gpp_hs_rpar_1m = read_nc(hs_rpar_1m, "mGPP") .* (3600 * 24 * 1e-6 * 12);
dgpp_1m = gpp_hs_epar_1m .- gpp_hs_rpar_1m;
dgpp_1s = zeros(360,180,4);
seasons = ["MAM", "JJA", "SON", "DJF"];
season_months = [[3,4,5], [6,7,8], [9,10,11], [12,1,2]];
for i in 1:4
    for ilon in 1:360, ilat in 1:180
        dgpp_1s[ilon,ilat,i] = nanmean(dgpp_1m[ilon,ilat,season_months[i]]);
    end;
    dgpp_1s[mask,i] .= 0;
end;

# 7. save the data
if !isfile("../output/8_global.nc")
    save_nc!("../output/8_global.nc", "SΔGPP", dgpp_1s, Dict{String,String}("about" => "Changes in GPP"));
    append_nc!("../output/8_global.nc", "ΔGPP", dgpp, Dict{String,String}("about" => "Changes in GPP"), ["lon", "lat"]);
    append_nc!("../output/8_global.nc", "RGPP", rgpp, Dict{String,String}("about" => "Relative changes in GPP relative to regular PAR defintion (%)"), ["lon", "lat"]);
    append_nc!("../output/8_global.nc", "mLAI", lai_1y, Dict{String,String}("about" => "1Y LAI from MODIS 8D 2019"), ["lon", "lat"]);
    append_nc!("../output/8_global.nc", "mCI", ci_1y, Dict{String,String}("about" => "1Y CI from MODIS"), ["lon", "lat"]);
    append_nc!("../output/8_global.nc", "mCHL", chl_1y, Dict{String,String}("about" => "1Y CHL from MODIS"), ["lon", "lat"]);
end;

# 8. display the fitting results
#=
lr_lai.R² = 0.8730182194194606
lr_ci.R² = 0.2334268237634859
lr_chl.R² = 0.2588099323728206
=#
lr_lai = linear_regress((lai_1y[dgpp .> 0][:],1), dgpp[dgpp .> 0][:]);
@show lr_lai.R²;
lr_ci = linear_regress((ci_1y[dgpp .> 0][:],1), dgpp[dgpp .> 0][:]);
@show lr_ci.R²;
lr_chl = linear_regress((chl_1y[dgpp .> 0][:],1), dgpp[dgpp .> 0][:]);
@show lr_chl.R²;

#=
lr_lai.R² = 0.8665152663007958
lr_ci.R² = 0.28729840856056954
lr_chl.R² = 0.1541600633669572
=#
lr_lai = linear_regress((lai_1y[rgpp .> 0][:],1), rgpp[rgpp .> 0][:]);
@show lr_lai.R²;
lr_ci = linear_regress((ci_1y[rgpp .> 0][:],1), rgpp[rgpp .> 0][:]);
@show lr_ci.R²;
lr_chl = linear_regress((chl_1y[rgpp .> 0][:],1), rgpp[rgpp .> 0][:]);
@show lr_chl.R²;
