#
# this script is meant to generate the surface albedo maps for the UV and PAR broadbands
#
using ProgressMeter: @showprogress

using Emerald.EmeraldLand.Namespace: BulkSPAC, LAND_2021_1NM, SPACConfiguration
using Emerald.EmeraldLand.SPAC: initialize_spac!, prescribe_traits!, soil_plant_air_continuum!
using Emerald.EmeraldMath.Stats: nanmax, nanmean
using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: leaf_chlorophyll_collection, query_collection
using GriddingMachine.Indexer: read_LUT
using NetcdfIO: append_nc!, save_nc!

FT = Float64;
config = SPACConfiguration(FT, dataset = LAND_2021_1NM);
config.ENABLE_REF = true;
config.ENABLE_SIF = false;
spac = BulkSPAC(config);
initialize_spac!(config, spac);
soil_plant_air_continuum!(config, spac, FT(0));


# 1. define the function to compute the albedos for UV and PAR broadbands
function uv_par_albedo_maps(chl_map::Matrix{FT}, lai_map::Matrix{FT}) where {FT<:AbstractFloat}
    (; Λ) = config.SPECTRA;
    mask_uv = (Λ .< 400);
    mask_par = (400 .<= Λ .<= 700);

    # compute the albedo for UV and PAR per chl+lai combo
    @inline uv_par_albedos(i) = (
        spac_tmp = deepcopy(spac);
        prescribe_traits!(config, spac_tmp; cab = chl_map[i], car = chl_map[i] / 7, lai = lai_map[i]);
        initialize_spac!(config, spac_tmp);
        soil_plant_air_continuum!(config, spac_tmp, FT(0));

        return nanmean(spac_tmp.canopy.sun_geometry.auxil.albedo[mask_uv]), nanmean(spac_tmp.canopy.sun_geometry.auxil.albedo[mask_par])
    );

    # create arrays of data, and fill them with data
    alb_uv = similar(chl_map) .* FT(NaN);
    alb_par = similar(chl_map) .* FT(NaN);
    @showprogress for i in eachindex(chl_map)
        if i%100 == 0
            @show i / length(chl_map) * 100;
        end;
        if chl_map[i]>0 && lai_map[i]>0
            alb_uv[i], alb_par[i] = uv_par_albedos(i);
        end;
    end;

    return alb_uv, alb_par
end;



# 2. create the ff maps based on global CHL and LAI maps
bias_nc = "../output/10_albedo.nc";
if !isfile(bias_nc)
    # generate the global mean CHL map
    chl_map = regrid(read_LUT(query_collection(leaf_chlorophyll_collection()))[1]);
    lai_map = regrid(read_LUT(query_collection("LAI_MODIS_2X_1M_2019_V1"))[1]);
    chl_ave = similar(chl_map[:,:,1]) .* FT(NaN);
    for i in axes(chl_map,1), j in axes(chl_map,2)
        chl_ave[i,j] = nanmean(chl_map[i,j,:]);
    end;
    lai_max = similar(lai_map[:,:,1]) .* FT(NaN);
    for i in axes(lai_map,1), j in axes(lai_map,2)
        lai_max[i,j] = nanmax(lai_map[i,j,:]);
    end;
    save_nc!(bias_nc, "CHL_MEAN", chl_ave, Dict{String,String}("about" => "Annually mean chlorophyll content"));
    append_nc!(bias_nc, "LAI_MAX", lai_max, Dict{String,String}("about" => "Maximum LAI"), ["lon", "lat"]);

    # generate the ff map for TOC
    uv_alb, par_alb = uv_par_albedo_maps(chl_ave, lai_max);
    append_nc!(bias_nc, "ALB_UV" , uv_alb , Dict{String,String}("about" => "Surface albedo for UV" ), ["lon", "lat"]);
    append_nc!(bias_nc, "ALB_PAR", par_alb, Dict{String,String}("about" => "Surface albedo for PAR"), ["lon", "lat"]);
end;
