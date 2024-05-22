#
# this script is meant to generate the sensitivity analysis of the PAR extension to
#     - Chl
#     - LAI
#     - CI
#     - Rad
#
using DataFrames: DataFrame

using Emerald.EmeraldIO.Text: save_csv!
using Emerald.EmeraldLand.Namespace: BulkSPAC, LAND_2021_1NM, ReferenceSpectra, SPACConfiguration
using Emerald.EmeraldLand.SPAC: GPP, initialize_spac!, prescribe_traits!, soil_plant_air_continuum!

FT = Float64;


# 1. create the config and spac
config_750 = SPACConfiguration(FT, dataset = LAND_2021_1NM, wl_par = [300,750], wl_par_700 = [300,700]);
config_750.ENABLE_REF = false;
config_750.ENABLE_SIF = false;
config_700 = SPACConfiguration(FT, dataset = LAND_2021_1NM, wl_par = [300,700], wl_par_700 = [300,700]);
config_700.ENABLE_REF = false;
config_700.ENABLE_SIF = false;
spac_750 = BulkSPAC(config_750);
spac_700 = BulkSPAC(config_700);
initialize_spac!(config_750, spac_750);
initialize_spac!(config_700, spac_700);
soil_plant_air_continuum!(config_750, spac_750, 1);
soil_plant_air_continuum!(config_700, spac_700, 1);
df = DataFrame();

# 2. create the sensitivity analysis to Chl
begin
    chls = collect(FT, 5:1:80);
    gpp_750 = FT[];
    gpp_700 = FT[];
    for i in eachindex(chls)
        spac_750_tmp = deepcopy(spac_750);
        spac_700_tmp = deepcopy(spac_700);
        prescribe_traits!(config_750, spac_750_tmp; cab = chls[i], car = chls[i] / 7);
        prescribe_traits!(config_700, spac_700_tmp; cab = chls[i], car = chls[i] / 7);
        initialize_spac!(config_750, spac_750_tmp);
        initialize_spac!(config_700, spac_700_tmp);
        for j in 1:360
            soil_plant_air_continuum!(config_750, spac_750_tmp, 10);
            soil_plant_air_continuum!(config_700, spac_700_tmp, 10);
        end;
        push!(gpp_750, GPP(spac_750_tmp));
        push!(gpp_700, GPP(spac_700_tmp));
    end;
    df[!, "CHL"] = chls;
    df[!, "GPP_750_CHL"] = gpp_750;
    df[!, "GPP_700_CHL"] = gpp_700;
end;

# 3. create the sensitivity analysis to LAI
begin
    lais = collect(FT, 0.5:0.1:8);
    gpp_750 = FT[];
    gpp_700 = FT[];
    for i in eachindex(lais)
        spac_750_tmp = deepcopy(spac_750);
        spac_700_tmp = deepcopy(spac_700);
        prescribe_traits!(config_750, spac_750_tmp; lai = lais[i]);
        prescribe_traits!(config_700, spac_700_tmp; lai = lais[i]);
        initialize_spac!(config_750, spac_750_tmp);
        initialize_spac!(config_700, spac_700_tmp);
        for j in 1:360
            soil_plant_air_continuum!(config_750, spac_750_tmp, 10);
            soil_plant_air_continuum!(config_700, spac_700_tmp, 10);
        end;
        push!(gpp_750, GPP(spac_750_tmp));
        push!(gpp_700, GPP(spac_700_tmp));
    end;
    df[!, "LAI"] = lais;
    df[!, "GPP_750_LAI"] = gpp_750;
    df[!, "GPP_700_LAI"] = gpp_700;
end;

# 4. create the sensitivity analysis to CI
begin
    clis = collect(FT, 0.25:0.01:1);
    gpp_750 = FT[];
    gpp_700 = FT[];
    for i in eachindex(clis)
        spac_750_tmp = deepcopy(spac_750);
        spac_700_tmp = deepcopy(spac_700);
        prescribe_traits!(config_750, spac_750_tmp; ci = clis[i]);
        prescribe_traits!(config_700, spac_700_tmp; ci = clis[i]);
        initialize_spac!(config_750, spac_750_tmp);
        initialize_spac!(config_700, spac_700_tmp);
        for j in 1:360
            soil_plant_air_continuum!(config_750, spac_750_tmp, 10);
            soil_plant_air_continuum!(config_700, spac_700_tmp, 10);
        end;
        push!(gpp_750, GPP(spac_750_tmp));
        push!(gpp_700, GPP(spac_700_tmp));
    end;
    df[!, "CI"] = clis;
    df[!, "GPP_750_CI"] = gpp_750;
    df[!, "GPP_700_CI"] = gpp_700;
end;

# 5. create the sensitivity analysis to SZA (and thus radiation)
begin
    szas = collect(FT, 0:1.14:85.51)[end:-1:1];
    gpp_750 = FT[];
    gpp_700 = FT[];
    rad_sw_bak = deepcopy(spac_750.meteo.rad_sw);
    for i in eachindex(szas)
        sza = szas[i];
        ratio = max(0, cosd(sza)) / cosd(48.317);
        spac_750_tmp = deepcopy(spac_750);
        spac_700_tmp = deepcopy(spac_700);
        spac_750_tmp.meteo.rad_sw.e_dir .= rad_sw_bak.e_dir .* ratio;
        spac_750_tmp.meteo.rad_sw.e_dif .= rad_sw_bak.e_dif .* ratio;
        spac_700_tmp.meteo.rad_sw.e_dir .= rad_sw_bak.e_dir .* ratio;
        spac_700_tmp.meteo.rad_sw.e_dif .= rad_sw_bak.e_dif .* ratio;
        spac_750_tmp.canopy.sun_geometry.state.sza = sza;
        spac_700_tmp.canopy.sun_geometry.state.sza = sza;
        for j in 1:360
            soil_plant_air_continuum!(config_750, spac_750_tmp, 10);
            soil_plant_air_continuum!(config_700, spac_700_tmp, 10);
        end;
        push!(gpp_750, GPP(spac_750_tmp));
        push!(gpp_700, GPP(spac_700_tmp));
    end;
    df[!, "SZA"] = szas;
    df[!, "GPP_750_SZA"] = gpp_750;
    df[!, "GPP_700_SZA"] = gpp_700;
end;

# 6. create the sensitivity analysis to SZA (and thus radiation)
begin
    szas = collect(FT, 0:1.14:85.51)[end:-1:1];
    gpp_750 = FT[];
    gpp_700 = FT[];
    rad_sw_bak = deepcopy(spac_750.meteo.rad_sw);
    for i in eachindex(szas)
        sza = szas[i];
        ratio = max(0, cosd(sza)) / cosd(48.317) * 2;
        spac_750_tmp = deepcopy(spac_750);
        spac_700_tmp = deepcopy(spac_700);
        spac_750_tmp.meteo.rad_sw.e_dir .= rad_sw_bak.e_dir .* ratio;
        spac_750_tmp.meteo.rad_sw.e_dif .= rad_sw_bak.e_dif .* ratio;
        spac_700_tmp.meteo.rad_sw.e_dir .= rad_sw_bak.e_dir .* ratio;
        spac_700_tmp.meteo.rad_sw.e_dif .= rad_sw_bak.e_dif .* ratio;
        spac_750_tmp.canopy.sun_geometry.state.sza = sza;
        spac_700_tmp.canopy.sun_geometry.state.sza = sza;
        for j in 1:360
            soil_plant_air_continuum!(config_750, spac_750_tmp, 10);
            soil_plant_air_continuum!(config_700, spac_700_tmp, 10);
        end;
        push!(gpp_750, GPP(spac_750_tmp));
        push!(gpp_700, GPP(spac_700_tmp));
    end;
    df[!, "SZA"] = szas;
    df[!, "GPP_750_SZA_2"] = gpp_750;
    df[!, "GPP_700_SZA_2"] = gpp_700;
end;

# 6. save the results
save_csv!("../output/7_sensitivity.csv", df);
