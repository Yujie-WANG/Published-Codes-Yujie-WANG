#
# this script is meant to generate the bias map as a proof of concept why we conduct the experiment at the global scale
#
using ProgressMeter: @showprogress

using Emerald.EmeraldLand.LeafOptics: leaf_spectra!
using Emerald.EmeraldLand.Namespace: BulkSPAC, LAND_2021_1NM, Leaf, SPACConfiguration
using Emerald.EmeraldLand.SPAC: initialize!, prescribe_traits!, soil_plant_air_continuum!
using Emerald.EmeraldMath.Stats: nanmax, nanmean
using Emerald.EmeraldPhysics.Optics: photon
using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: leaf_chlorophyll_collection, query_collection
using GriddingMachine.Indexer: read_LUT
using NetcdfIO: append_nc!, save_nc!

FT = Float64;
config = SPACConfiguration{FT}(DATASET = LAND_2021_1NM, ENABLE_REF = false, ENABLE_SIF = false);
leaf = Leaf(config);


# 1. define the faparfppar map function
function f_apar_f_ppar(chl_map::Matrix{FT}, lai_map::Matrix{FT}) where {FT<:AbstractFloat}
    (; IΛ_PAR, ΔΛ_PAR, Λ) = config.SPECTRA;

    # compute the ratios for CliMA Land PRO
    @inline ratios(chl::Number, rad::Vector) = (
        leaf.bio.state.cab = chl;
        leaf.bio.state.car = chl / 7;
        leaf_spectra!(config, leaf.bio, FT(5));

        photon_all = photon.(Λ, rad) .* 1000;
        par = photon_all[IΛ_PAR]' * ΔΛ_PAR;

        photon_ppar = leaf.bio.auxil.f_ppar .* leaf.bio.auxil.α_leaf .* photon_all;
        ppar = photon_ppar[IΛ_PAR]' * ΔΛ_PAR;

        return ppar / par
    );

    # create arrays of data, and fill them with data
    map_ff_toc = similar(chl_map) .* FT(NaN);
    map_ff_boc = similar(chl_map) .* FT(NaN);
    @showprogress for i in eachindex(chl_map)
        if chl_map[i]>0 && lai_map[i]>0
            spac = BulkSPAC(config);
            prescribe_traits!(config, spac; cab = chl_map[i], car = chl_map[i] / 7, lai = lai_map[i]);
            initialize!(config, spac);
            soil_plant_air_continuum!(config, spac, FT(0));
            rad = spac.canopy.sun_geometry.auxil.e_dirꜜ[:,1] .+ spac.canopy.sun_geometry.auxil.e_difꜜ[:,1] .+ spac.canopy.sun_geometry.auxil.e_difꜛ[:,1];
            map_ff_toc[i] = ratios(chl_map[i], rad);

            rad = spac.canopy.sun_geometry.auxil.e_dirꜜ[:,end] .+ spac.canopy.sun_geometry.auxil.e_difꜜ[:,end] .+ spac.canopy.sun_geometry.auxil.e_difꜛ[:,end];
            map_ff_boc[i] = ratios(chl_map[i], rad);
        end;
    end;

    return map_ff_toc, map_ff_boc
end;

# 2. create the ff maps based on global CHL and LAI maps
bias_nc = "../output/3_bias.nc";
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
    spac = BulkSPAC(config);
    prescribe_traits!(config, spac; lai = 3);
    initialize!(config, spac);
    soil_plant_air_continuum!(config, spac, FT(0));
    rad = spac.canopy.sun_geometry.auxil.e_dirꜜ[:,1] .+ spac.canopy.sun_geometry.auxil.e_difꜜ[:,1] .+ spac.canopy.sun_geometry.auxil.e_difꜛ[:,1];
    ff_toc, ff_boc = f_apar_f_ppar(chl_ave, lai_max);
    append_nc!(bias_nc, "FF_TOC", ff_toc, Dict{String,String}("about" => "F_APAR * F_PPAR based on global mean chlorophyll content at TOC"), ["lon", "lat"]);
    append_nc!(bias_nc, "FF_BOC", ff_boc, Dict{String,String}("about" => "F_APAR * F_PPAR based on global mean chlorophyll content at BOC"), ["lon", "lat"]);
end;
