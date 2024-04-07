#
# this script is used to create a sensitivity analysis of the leaf's absorption features, PAR, and A
#
using DataFrames: DataFrame

using Emerald.EmeraldIO.Text: save_csv!
using Emerald.EmeraldLand.Namespace: BulkSPAC, LAND_2021_1NM, Leaf, SPACConfiguration
using Emerald.EmeraldLand.SPAC: initialize!, prescribe_traits!, soil_plant_air_continuum!
using Emerald.EmeraldPhysics.Optics: photon

FT = Float64;


# 1. define function to compute PAR, APAR, and PPAR from Leaf and ShortwaveRadiation
function leaf_fapar_fppar(config::SPACConfiguration{FT}, leaf::Leaf{FT}, rad::Vector{FT}) where {FT}
    (; IΛ_PAR, ΔΛ_PAR, Λ) = config.SPECTRA;

    # convert the radiation to photon flux
    photon_all = photon.(Λ, rad) .* 1000;

    # compute the absorbed PAR by the leaf
    photon_leaf = leaf.bio.auxil.α_leaf .* photon_all;

    # compute the PAR that goes to photosystems
    photon_ppar = leaf.bio.auxil.f_ppar .* photon_leaf;

    # compute the sum of par, apar, and ppar
    par = photon_all[IΛ_PAR]' * ΔΛ_PAR;
    apar = photon_leaf[IΛ_PAR]' * ΔΛ_PAR;
    ppar = photon_ppar[IΛ_PAR]' * ΔΛ_PAR;

    return apar / par, ppar / apar
end;

# 2. create the configuration and spac structs for the simulation at the default radiation
config = SPACConfiguration{FT}(DATASET = LAND_2021_1NM, ENABLE_SIF = false);
spac = BulkSPAC(config);
initialize!(config, spac);
soil_plant_air_continuum!(config, spac, FT(0));
df = DataFrame(FAPAR = zeros(FT,length(spac.plant.leaves)),
               FAPAR_CAR = zeros(FT,length(spac.plant.leaves)),
               FAPAR_CHL = zeros(FT,length(spac.plant.leaves)),
               FAPAR_CI = zeros(FT,length(spac.plant.leaves)),
               FAPAR_LAI = zeros(FT,length(spac.plant.leaves)),
               FAPAR_LMA = zeros(FT,length(spac.plant.leaves)),
               FPPAR = zeros(FT,length(spac.plant.leaves)),
               FPPAR_CAR = zeros(FT,length(spac.plant.leaves)),
               FPPAR_CHL = zeros(FT,length(spac.plant.leaves)),
               FPPAR_CI = zeros(FT,length(spac.plant.leaves)),
               FPPAR_LAI = zeros(FT,length(spac.plant.leaves)),
               FPPAR_LMA = zeros(FT,length(spac.plant.leaves)));
for ilf in eachindex(spac.plant.leaves)[end:-1:1]
    irt = length(spac.plant.leaves) - ilf + 1;
    rad = spac.canopy.sun_geometry.auxil.e_dirꜜ[:,irt] .+ spac.canopy.sun_geometry.auxil.e_difꜜ[:,irt] .+ spac.canopy.sun_geometry.auxil.e_difꜛ[:,irt];
    cleaf = spac.plant.leaves[ilf];
    df.FAPAR[irt], df.FPPAR[irt] = leaf_fapar_fppar(config, cleaf, rad);
end;

# 3. do the simulation for the case with higher CHL
spac_chl = deepcopy(spac);
prescribe_traits!(config, spac_chl; cab = 80);
soil_plant_air_continuum!(config, spac_chl, FT(0));
for ilf in eachindex(spac_chl.plant.leaves)[end:-1:1]
    irt = length(spac_chl.plant.leaves) - ilf + 1;
    rad = spac_chl.canopy.sun_geometry.auxil.e_dirꜜ[:,irt] .+ spac_chl.canopy.sun_geometry.auxil.e_difꜜ[:,irt] .+ spac_chl.canopy.sun_geometry.auxil.e_difꜛ[:,irt];
    cleaf = spac_chl.plant.leaves[ilf];
    df.FAPAR_CHL[irt], df.FPPAR_CHL[irt] = leaf_fapar_fppar(config, cleaf, rad);
end;

# 4. do the simulation for the case with higher CAR
spac_car = deepcopy(spac);
prescribe_traits!(config, spac_car; car = 80 / 7);
soil_plant_air_continuum!(config, spac_car, FT(0));
for ilf in eachindex(spac_car.plant.leaves)[end:-1:1]
    irt = length(spac_car.plant.leaves) - ilf + 1;
    rad = spac_car.canopy.sun_geometry.auxil.e_dirꜜ[:,irt] .+ spac_car.canopy.sun_geometry.auxil.e_difꜜ[:,irt] .+ spac_car.canopy.sun_geometry.auxil.e_difꜛ[:,irt];
    cleaf = spac_car.plant.leaves[ilf];
    df.FAPAR_CAR[irt], df.FPPAR_CAR[irt] = leaf_fapar_fppar(config, cleaf, rad);
end;

# 5. do the simulation for the case with higher LMA
spac_lma = deepcopy(spac);
for cleaf in spac_lma.plant.leaves
    cleaf.bio.state.lma = 0.024;
end;
initialize!(config, spac_lma);
soil_plant_air_continuum!(config, spac_lma, FT(0));
for ilf in eachindex(spac_lma.plant.leaves)[end:-1:1]
    irt = length(spac_lma.plant.leaves) - ilf + 1;
    rad = spac_lma.canopy.sun_geometry.auxil.e_dirꜜ[:,irt] .+ spac_lma.canopy.sun_geometry.auxil.e_difꜜ[:,irt] .+ spac_lma.canopy.sun_geometry.auxil.e_difꜛ[:,irt];
    cleaf = spac_lma.plant.leaves[ilf];
    df.FAPAR_LMA[irt], df.FPPAR_LMA[irt] = leaf_fapar_fppar(config, cleaf, rad);
end;

# 6. do the simulation for the case with higher LAI
spac_lai = deepcopy(spac);
prescribe_traits!(config, spac_lai; lai = 6);
soil_plant_air_continuum!(config, spac_lai, FT(0));
for ilf in eachindex(spac_lai.plant.leaves)[end:-1:1]
    irt = length(spac_lai.plant.leaves) - ilf + 1;
    rad = spac_lai.canopy.sun_geometry.auxil.e_dirꜜ[:,irt] .+ spac_lai.canopy.sun_geometry.auxil.e_difꜜ[:,irt] .+ spac_lai.canopy.sun_geometry.auxil.e_difꜛ[:,irt];
    cleaf = spac_lai.plant.leaves[ilf];
    df.FAPAR_LAI[irt], df.FPPAR_LAI[irt] = leaf_fapar_fppar(config, cleaf, rad);
end;

# 7. do the simulation for the case with more heterogeneous horizontal structure (lower ci)
spac_ci = deepcopy(spac);
prescribe_traits!(config, spac_ci; ci = 0.6);
soil_plant_air_continuum!(config, spac_ci, FT(0));
for ilf in eachindex(spac_ci.plant.leaves)[end:-1:1]
    irt = length(spac_ci.plant.leaves) - ilf + 1;
    rad = spac_ci.canopy.sun_geometry.auxil.e_dirꜜ[:,irt] .+ spac_ci.canopy.sun_geometry.auxil.e_difꜜ[:,irt] .+ spac_ci.canopy.sun_geometry.auxil.e_difꜛ[:,irt];
    cleaf = spac_ci.plant.leaves[ilf];
    df.FAPAR_CI[irt], df.FPPAR_CI[irt] = leaf_fapar_fppar(config, cleaf, rad);
end;

# 8. save the results
save_csv!("../output/2_sensitivity.csv", df);
