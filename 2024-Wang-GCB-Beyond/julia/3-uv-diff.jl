#
# This script is meant to generate the model output to illustrate the impact of adding UV radiation to the model
#
using DataFrames: DataFrame

using Emerald.EmeraldIO.Text: save_csv!
using Emerald.EmeraldLand.Namespace: BulkSPAC, C3CytoState, C3CytoTrait, LAND_2021_1NM, SPACConfiguration
using Emerald.EmeraldLand.SPAC: initialize_spac!, soil_plant_air_continuum!
using Emerald.EmeraldMath.Stats: nanmean

FT = Float64;


# 1. create the configuration and spac structs for the simulation
config = SPACConfiguration(FT, dataset = LAND_2021_1NM);
spac_uv = BulkSPAC(config);
spac_uv_cyto = deepcopy(spac_uv);
for l in spac_uv_cyto.plant.leaves
    l.photosystem.trait = C3CytoTrait{FT}();
    l.photosystem.state = C3CytoState{FT}();
end;
spac_nouv = deepcopy(spac_uv);
spac_nouv_cyto = deepcopy(spac_uv_cyto);

# 2. mask out the UV radiation for the spac_nouv
mask_uv = config.SPECTRA.Λ .< 400;
spac_nouv.meteo.rad_sw.e_dir[mask_uv] .= 0;
spac_nouv.meteo.rad_sw.e_dif[mask_uv] .= 0;
spac_nouv_cyto.meteo.rad_sw.e_dir[mask_uv] .= 0;
spac_nouv_cyto.meteo.rad_sw.e_dif[mask_uv] .= 0;

initialize_spac!(config, spac_uv);
initialize_spac!(config, spac_uv_cyto);
initialize_spac!(config, spac_nouv);
initialize_spac!(config, spac_nouv_cyto);
soil_plant_air_continuum!(config, spac_uv, 1);
soil_plant_air_continuum!(config, spac_uv_cyto, 1);
soil_plant_air_continuum!(config, spac_nouv, 1);
soil_plant_air_continuum!(config, spac_nouv_cyto, 1);

# 3. test the case with/without UV radiation
for i in 1:360
    soil_plant_air_continuum!(config, spac_uv, 10);
    soil_plant_air_continuum!(config, spac_uv_cyto, 10);
    soil_plant_air_continuum!(config, spac_nouv, 10);
    soil_plant_air_continuum!(config, spac_nouv_cyto, 10);
end;
es_uv = FT[spac_uv.canopy.sun_geometry.auxil.r_net_sw_leaf .+ spac_uv.canopy.sun_geometry.auxil.r_net_sw_stem; spac_uv.soil_bulk.auxil.r_net_sw];
es_uv_cyto = FT[spac_uv_cyto.canopy.sun_geometry.auxil.r_net_sw_leaf .+ spac_uv_cyto.canopy.sun_geometry.auxil.r_net_sw_stem; spac_uv_cyto.soil_bulk.auxil.r_net_sw];
es_nouv = FT[spac_nouv.canopy.sun_geometry.auxil.r_net_sw_leaf .+ spac_nouv.canopy.sun_geometry.auxil.r_net_sw_stem; spac_nouv.soil_bulk.auxil.r_net_sw];
es_nouv_cyto = FT[spac_nouv_cyto.canopy.sun_geometry.auxil.r_net_sw_leaf .+ spac_nouv_cyto.canopy.sun_geometry.auxil.r_net_sw_stem; spac_nouv_cyto.soil_bulk.auxil.r_net_sw];
ts_uv = FT[[l.energy.s_aux.t for l in spac_uv.plant.leaves[end:-1:1]]; spac_uv.soils[1].s_aux.t];
ts_uv_cyto = FT[[l.energy.s_aux.t for l in spac_uv_cyto.plant.leaves[end:-1:1]]; spac_uv_cyto.soils[1].s_aux.t];
ts_nouv = FT[[l.energy.s_aux.t for l in spac_nouv.plant.leaves[end:-1:1]]; spac_nouv.soils[1].s_aux.t];
ts_nouv_cyto = FT[[l.energy.s_aux.t for l in spac_nouv_cyto.plant.leaves[end:-1:1]]; spac_nouv_cyto.soils[1].s_aux.t];
as_uv = FT[];
as_uv_cyto = FT[];
as_nouv = FT[];
as_nouv_cyto = FT[];
for i in eachindex(spac_uv.plant.leaves)[end:-1:1]
    leaf_uv = spac_uv.plant.leaves[i];
    leaf_uv_cyto = spac_uv_cyto.plant.leaves[i];
    leaf_nouv = spac_nouv.plant.leaves[i];
    leaf_nouv_cyto = spac_nouv_cyto.plant.leaves[i];
    f_sunlit = spac_uv.canopy.sun_geometry.s_aux.p_sunlit[length(spac_uv.plant.leaves)+1-i];
    a_uv = nanmean(leaf_uv.flux.auxil.a_g_sunlit) * f_sunlit + leaf_uv.flux.auxil.a_g_shaded * (1 - f_sunlit);
    a_uv_cyto = nanmean(leaf_uv_cyto.flux.auxil.a_g_sunlit) * f_sunlit + leaf_uv_cyto.flux.auxil.a_g_shaded * (1 - f_sunlit);
    a_nouv = nanmean(leaf_nouv.flux.auxil.a_g_sunlit) * f_sunlit + leaf_nouv.flux.auxil.a_g_shaded * (1 - f_sunlit);
    a_nouv_cyto = nanmean(leaf_nouv_cyto.flux.auxil.a_g_sunlit) * f_sunlit + leaf_nouv_cyto.flux.auxil.a_g_shaded * (1 - f_sunlit);
    push!(as_uv, a_uv);
    push!(as_uv_cyto, a_uv_cyto);
    push!(as_nouv, a_nouv);
    push!(as_nouv_cyto, a_nouv_cyto);
end;
push!(as_uv, NaN);
push!(as_uv_cyto, NaN);
push!(as_nouv, NaN);
push!(as_nouv_cyto, NaN);

# 4. save the results
save_csv!("../output/3_uv_diff.csv",
          DataFrame(ES_UV = es_uv,
                    ES_UV_CYTO = es_uv_cyto,
                    ES_NOUV = es_nouv,
                    ES_NOUV_CYTO = es_nouv_cyto,
                    TS_UV = ts_uv,
                    TS_UV_CYTO = ts_uv_cyto,
                    TS_NOUV = ts_nouv,
                    TS_NOUV_CYTO = ts_nouv_cyto,
                    AS_UV = as_uv,
                    AS_UV_CYTO = as_uv_cyto,
                    AS_NOUV = as_nouv,
                    AS_NOUV_CYTO = as_nouv_cyto));
