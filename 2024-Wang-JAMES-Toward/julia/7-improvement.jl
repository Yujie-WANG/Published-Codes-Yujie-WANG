#
# this script is meant to compare whether the use of the new empirical function improves the GPP estimation
#
using DataFrames: DataFrame
using ProgressMeter: @showprogress

using Emerald.EmeraldIO.Text: save_csv!
using Emerald.EmeraldLand.Namespace: BulkSPAC, LAND_2021_1NM, SPACConfiguration
using Emerald.EmeraldLand.SPAC: GPP, TROPOMI_SIF740, initialize!, prescribe_traits!, soil_plant_air_continuum!
using Emerald.EmeraldMath.Stats: nanmean
using Emerald.EmeraldPhysics.EarthGeometry: solar_zenith_angle

FT = Float64;


# 1. create spac with default setup
config = SPACConfiguration{FT}(DATASET = LAND_2021_1NM);
spac_def_1 = BulkSPAC(config);
spac_def_2 = BulkSPAC(config);
prescribe_traits!(config, spac_def_1; lai = 4, cab = 60, car = 60 / 7, vcmax = 60);
prescribe_traits!(config, spac_def_2; lai = 4, cab = 60, car = 60 / 7, vcmax = 30);
initialize!(config, spac_def_1);
initialize!(config, spac_def_2);
soil_plant_air_continuum!(config, spac_def_1, FT(0));
soil_plant_air_continuum!(config, spac_def_2, FT(0));

# 2. create spac with CLM5 set up for BET tropical forest
spac_clm_1 = deepcopy(spac_def_1);
spac_clm_2 = deepcopy(spac_def_2);
for leaf in spac_clm_1.plant.leaves
    leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_PAR_700] .= 0.10;
    leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_PAR_700] .= 0.05;
    leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_NIR] .= 0.45;
    leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_NIR] .= 0.25;
    leaf.bio.auxil.α_leaf .= 1 .- leaf.bio.auxil.ρ_leaf .- leaf.bio.auxil.τ_leaf;
    leaf.bio.auxil.f_ppar[config.SPECTRA.IΛ_PAR_700] .= 1;
    leaf.bio.auxil.f_sife[config.SPECTRA.IΛ_PAR_700] .= 1;
end;
for leaf in spac_clm_2.plant.leaves
    leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_PAR_700] .= 0.07;
    leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_PAR_700] .= 0.05;
    leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_NIR] .= 0.35;
    leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_NIR] .= 0.10;
    leaf.bio.auxil.α_leaf .= 1 .- leaf.bio.auxil.ρ_leaf .- leaf.bio.auxil.τ_leaf;
    leaf.bio.auxil.f_ppar[config.SPECTRA.IΛ_PAR_700] .= 1;
    leaf.bio.auxil.f_sife[config.SPECTRA.IΛ_PAR_700] .= 1;
end;

# 3. create spac with empirical ρ, τ, and f_ppar corrections
spac_emp_1 = deepcopy(spac_def_1);
spac_emp_2 = deepcopy(spac_def_2);
for ilf in eachindex(spac_emp_1.plant.leaves)
    leaf = spac_emp_1.plant.leaves[ilf];
    f_ppar = -0.0374643 * ilf / length(spac_emp_2.plant.leaves) + -0.0009303 * spac_emp_2.canopy.structure.state.lai + 0.193018 * log(log(leaf.bio.state.cab)) + 0.529997;
    ρ_par = nanmean(leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_PAR_700]);
    τ_par = nanmean(leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_PAR_700]);
    ρ_nir = nanmean(leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_NIR]);
    τ_nir = nanmean(leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_NIR]);
    leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_PAR_700] .= ρ_par;
    leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_PAR_700] .= τ_par;
    leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_NIR] .= ρ_nir;
    leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_NIR] .= τ_nir;
    leaf.bio.auxil.α_leaf .= 1 .- leaf.bio.auxil.ρ_leaf .- leaf.bio.auxil.τ_leaf;
    leaf.bio.auxil.f_ppar[config.SPECTRA.IΛ_PAR_700] .= f_ppar;
    leaf.bio.auxil.f_sife[config.SPECTRA.IΛ_PAR_700] .= f_ppar;
end;
for ilf in eachindex(spac_emp_2.plant.leaves)
    leaf = spac_emp_2.plant.leaves[ilf];
    f_ppar = -0.0374643 * ilf / length(spac_emp_2.plant.leaves) + -0.0009303 * spac_emp_2.canopy.structure.state.lai + 0.193018 * log(log(leaf.bio.state.cab)) + 0.529997;
    ρ_par = nanmean(leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_PAR_700]);
    τ_par = nanmean(leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_PAR_700]);
    ρ_nir = nanmean(leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_NIR]);
    τ_nir = nanmean(leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_NIR]);
    leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_PAR_700] .= ρ_par;
    leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_PAR_700] .= τ_par;
    leaf.bio.auxil.ρ_leaf[config.SPECTRA.IΛ_NIR] .= ρ_nir;
    leaf.bio.auxil.τ_leaf[config.SPECTRA.IΛ_NIR] .= τ_nir;
    leaf.bio.auxil.α_leaf .= 1 .- leaf.bio.auxil.ρ_leaf .- leaf.bio.auxil.τ_leaf;
    leaf.bio.auxil.f_ppar[config.SPECTRA.IΛ_PAR_700] .= f_ppar;
    leaf.bio.auxil.f_sife[config.SPECTRA.IΛ_PAR_700] .= f_ppar;
end;

# 4. run a diurnal cycle with each spac
#    TODO: when sza is 90, there could be a problem with the soil_plant_air_continuum! function (FIX IT)
df = DataFrame(Hour = Float64[],
               GPP_DEF_1 = Float64[],
               GPP_CLM_1 = Float64[],
               GPP_EMP_1 = Float64[],
               GPP_DEF_2 = Float64[],
               GPP_CLM_2 = Float64[],
               GPP_EMP_2 = Float64[],
               SIF_DEF_1 = Float64[],
               SIF_CLM_1 = Float64[],
               SIF_EMP_1 = Float64[],
               SIF_DEF_2 = Float64[],
               SIF_CLM_2 = Float64[],
               SIF_EMP_2 = Float64[]);
rad_bak = deepcopy(spac_def_1.meteo.rad_sw);
@showprogress for hour in collect(0.001:0.1:24)
    sza_1 = solar_zenith_angle(0.0, 180.0, FT(hour), FT(0));
    spac_def_1.canopy.sun_geometry.state.sza = sza_1;
    spac_clm_1.canopy.sun_geometry.state.sza = sza_1;
    spac_emp_1.canopy.sun_geometry.state.sza = sza_1;

    sza_2 = solar_zenith_angle(45.0, 180.0, FT(hour), FT(0));
    spac_def_2.canopy.sun_geometry.state.sza = sza_2;
    spac_clm_2.canopy.sun_geometry.state.sza = sza_2;
    spac_emp_2.canopy.sun_geometry.state.sza = sza_2;

    ratio_1 = max(0, cosd(sza_1)) / cosd(48.317);
    spac_def_1.meteo.rad_sw.e_dir .= rad_bak.e_dir .* ratio_1;
    spac_def_1.meteo.rad_sw.e_dif .= rad_bak.e_dif .* ratio_1;
    spac_clm_1.meteo.rad_sw.e_dir .= rad_bak.e_dir .* ratio_1;
    spac_clm_1.meteo.rad_sw.e_dif .= rad_bak.e_dif .* ratio_1;
    spac_emp_1.meteo.rad_sw.e_dir .= rad_bak.e_dir .* ratio_1;
    spac_emp_1.meteo.rad_sw.e_dif .= rad_bak.e_dif .* ratio_1;

    ratio_2 = max(0, cosd(sza_2)) / cosd(48.317);
    spac_def_2.meteo.rad_sw.e_dir .= rad_bak.e_dir .* ratio_2;
    spac_def_2.meteo.rad_sw.e_dif .= rad_bak.e_dif .* ratio_2;
    spac_clm_2.meteo.rad_sw.e_dir .= rad_bak.e_dir .* ratio_2;
    spac_clm_2.meteo.rad_sw.e_dif .= rad_bak.e_dif .* ratio_2;
    spac_emp_2.meteo.rad_sw.e_dir .= rad_bak.e_dir .* ratio_2;
    spac_emp_2.meteo.rad_sw.e_dif .= rad_bak.e_dif .* ratio_2;

    for i in 1:6
        soil_plant_air_continuum!(config, spac_def_1, FT(60));
        soil_plant_air_continuum!(config, spac_clm_1, FT(60));
        soil_plant_air_continuum!(config, spac_emp_1, FT(60));
        soil_plant_air_continuum!(config, spac_def_2, FT(60));
        soil_plant_air_continuum!(config, spac_clm_2, FT(60));
        soil_plant_air_continuum!(config, spac_emp_2, FT(60));
    end;

    push!(df, (hour,
               GPP(spac_def_1),
               GPP(spac_clm_1),
               GPP(spac_emp_1),
               GPP(spac_def_2),
               GPP(spac_clm_2),
               GPP(spac_emp_2),
               TROPOMI_SIF740(config, spac_def_1),
               TROPOMI_SIF740(config, spac_clm_1),
               TROPOMI_SIF740(config, spac_emp_1),
               TROPOMI_SIF740(config, spac_def_2),
               TROPOMI_SIF740(config, spac_clm_2),
               TROPOMI_SIF740(config, spac_emp_2)));
end;

# 5. save the results to a csv file
save_csv!("../output/7_improvement.csv", df);
