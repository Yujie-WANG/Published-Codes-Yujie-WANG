#
# this script is used to output the model simulation for the new far-red radiation feature
#
using DataFrames: DataFrame
using ProgressMeter: @showprogress
using Statistics: mean

using Emerald.EmeraldIO.Text: save_csv!
using Emerald.EmeraldLand.Namespace: BulkSPAC, C3CytoState, C3CytoTrait, LAND_2021_1NM, ReferenceSpectra, SPACConfiguration
using Emerald.EmeraldLand.SPAC: GPP, initialize_spac!, soil_plant_air_continuum!

FT = Float64;


# 1. sensitivity to change in radiation to show that the model is working
begin
    # 1.1 create the configuration and spac structs for the simulation
    config = SPACConfiguration(FT, dataset = LAND_2021_1NM);
    config.ENABLE_REF = false;
    config.ENABLE_SIF = false;
    spac = BulkSPAC(config);
    spac_cyto = deepcopy(spac);
    for l in spac_cyto.plant.leaves
        l.photosystem.trait = C3CytoTrait{FT}();
        l.photosystem.state = C3CytoState{FT}();
    end;
    initialize_spac!(config, spac);
    initialize_spac!(config, spac_cyto);
    soil_plant_air_continuum!(config, spac, FT(100));
    soil_plant_air_continuum!(config, spac_cyto, FT(100));

    # 1.2 set up the masks to change the radiation
    mask_farr = 700 .< config.SPECTRA.Λ .<= 750;
    mask_rpar = 400 .<= config.SPECTRA.Λ .<= 700;
    mask_uvio = config.SPECTRA.Λ .< 400;

    # 1.3 change the radiation and run the simulation
    df = DataFrame();
    for label in ["Ratio", "RPAR_UV", "RPAR_UV_FR", "FR", "2%RPAR_FR", "RPAR_UV_CYTO", "RPAR_UV_FR_CYTO", "FR_CYTO", "2%RPAR_FR_CYTO"]
        df[!, label] = Float64[];
    end;
    @showprogress for ratio in [0.01:0.01:0.1; 0.15:0.05:1]
        # spac with rPAR and UV
        spac_1 = deepcopy(spac);
        spac_1.meteo.rad_sw.e_dir .*= ratio;
        spac_1.meteo.rad_sw.e_dif .*= ratio;
        spac_1.meteo.rad_sw.e_dir[mask_farr] .= 0;
        spac_1.meteo.rad_sw.e_dif[mask_farr] .= 0;
        spac_1_cyto = deepcopy(spac_cyto);
        spac_1_cyto.meteo.rad_sw.e_dir .*= ratio;
        spac_1_cyto.meteo.rad_sw.e_dif .*= ratio;
        spac_1_cyto.meteo.rad_sw.e_dir[mask_farr] .= 0;
        spac_1_cyto.meteo.rad_sw.e_dif[mask_farr] .= 0;
        for j in 1:60
            soil_plant_air_continuum!(config, spac_1, 60);
            soil_plant_air_continuum!(config, spac_1_cyto, 60);
        end;

        # spac with rPAR, UV, and FR
        spac_3 = deepcopy(spac);
        spac_3.meteo.rad_sw.e_dir .*= ratio;
        spac_3.meteo.rad_sw.e_dif .*= ratio;
        spac_3_cyto = deepcopy(spac_cyto);
        spac_3_cyto.meteo.rad_sw.e_dir .*= ratio;
        spac_3_cyto.meteo.rad_sw.e_dif .*= ratio;
        for j in 1:60
            soil_plant_air_continuum!(config, spac_3, 60);
            soil_plant_air_continuum!(config, spac_3_cyto, 60);
        end;

        # spac with FR only
        spac_4 = deepcopy(spac);
        spac_4.meteo.rad_sw.e_dir .*= ratio;
        spac_4.meteo.rad_sw.e_dif .*= ratio;
        spac_4.meteo.rad_sw.e_dir[mask_rpar] .= 0;
        spac_4.meteo.rad_sw.e_dif[mask_rpar] .= 0;
        spac_4.meteo.rad_sw.e_dir[mask_uvio] .= 0;
        spac_4.meteo.rad_sw.e_dif[mask_uvio] .= 0;
        spac_4_cyto = deepcopy(spac_cyto);
        spac_4_cyto.meteo.rad_sw.e_dir .*= ratio;
        spac_4_cyto.meteo.rad_sw.e_dif .*= ratio;
        spac_4_cyto.meteo.rad_sw.e_dir[mask_rpar] .= 0;
        spac_4_cyto.meteo.rad_sw.e_dif[mask_rpar] .= 0;
        spac_4_cyto.meteo.rad_sw.e_dir[mask_uvio] .= 0;
        spac_4_cyto.meteo.rad_sw.e_dif[mask_uvio] .= 0;
        for j in 1:60
            soil_plant_air_continuum!(config, spac_4, 60);
            soil_plant_air_continuum!(config, spac_4_cyto, 60);
        end;

        # spac with 2% of rPAR and ratio of FR
        spac_5 = deepcopy(spac);
        spac_5.meteo.rad_sw.e_dir .*= ratio;
        spac_5.meteo.rad_sw.e_dif .*= ratio;
        spac_5.meteo.rad_sw.e_dir[mask_uvio .|| mask_rpar] .= spac.meteo.rad_sw.e_dir[mask_uvio .|| mask_rpar] .* 0.02;
        spac_5.meteo.rad_sw.e_dif[mask_uvio .|| mask_rpar] .= spac.meteo.rad_sw.e_dif[mask_uvio .|| mask_rpar] .* 0.02
        spac_5_cyto = deepcopy(spac_cyto);
        spac_5_cyto.meteo.rad_sw.e_dir .*= ratio;
        spac_5_cyto.meteo.rad_sw.e_dif .*= ratio;
        spac_5_cyto.meteo.rad_sw.e_dir[mask_uvio .|| mask_rpar] .= spac.meteo.rad_sw.e_dir[mask_uvio .|| mask_rpar] .* 0.02;
        spac_5_cyto.meteo.rad_sw.e_dif[mask_uvio .|| mask_rpar] .= spac.meteo.rad_sw.e_dif[mask_uvio .|| mask_rpar] .* 0.02
        for j in 1:60
            soil_plant_air_continuum!(config, spac_5, 60);
            soil_plant_air_continuum!(config, spac_5_cyto, 60);
        end;

        push!(df, [ratio, GPP(spac_1), GPP(spac_3), GPP(spac_4), GPP(spac_5), GPP(spac_1_cyto), GPP(spac_3_cyto), GPP(spac_4_cyto), GPP(spac_5_cyto)]);
    end;

    # 1.4 save the data
    save_csv!("../output/5_fr_diff.csv", df);
end;

# 2. get the profiles within the canopy
begin
    # 2.1 run the model at two contrasting configurations
    config_750 = SPACConfiguration(FT, dataset = LAND_2021_1NM, wl_par = [300,750], wl_par_700 = [300,700]);
    config_750.ENABLE_REF = false;
    config_750.ENABLE_SIF = false;
    config_700 = SPACConfiguration(FT, dataset = LAND_2021_1NM, wl_par = [300,700], wl_par_700 = [300,700]);
    config_700.ENABLE_REF = false;
    config_700.ENABLE_SIF = false;
    spac_750 = BulkSPAC(config_750);
    spac_700 = BulkSPAC(config_700);
    spac_750_cyto = deepcopy(spac_750);
    spac_700_cyto = deepcopy(spac_700);
    for l in spac_750_cyto.plant.leaves
        l.photosystem.trait = C3CytoTrait{FT}();
        l.photosystem.state = C3CytoState{FT}();
    end;
    for l in spac_700_cyto.plant.leaves
        l.photosystem.trait = C3CytoTrait{FT}();
        l.photosystem.state = C3CytoState{FT}();
    end;
    initialize_spac!(config_750, spac_750);
    initialize_spac!(config_700, spac_700);
    initialize_spac!(config_750, spac_750_cyto);
    initialize_spac!(config_700, spac_700_cyto);
    for j in 1:60
        soil_plant_air_continuum!(config_750, spac_750, 60);
        soil_plant_air_continuum!(config_700, spac_700, 60);
        soil_plant_air_continuum!(config_750, spac_750_cyto, 60);
        soil_plant_air_continuum!(config_700, spac_700_cyto, 60);
    end;

    # 2.2 get the profiles for UV, PAR, and FR, Ag, and gs
    df = DataFrame();
    for label in ["UV", "PAR", "FR", "GREEN", "AG_750", "AG_700", "GS_750", "GS_700", "UV_CYTO", "PAR_CYTO", "FR_CYTO", "GREEN_CYTO", "AG_750_CYTO", "AG_700_CYTO", "GS_750_CYTO", "GS_700_CYTO"]
        df[!, label] = Float64[];
    end;
    mask_grin = 500 .<= config_750.SPECTRA.Λ .<= 570;
    for ilf in eachindex(spac_750.plant.leaves)[end:-1:1]
        irt = length(spac_750.plant.leaves) - ilf + 1;
        rad = spac_750.canopy.sun_geometry.auxil.e_dirꜜ[:,irt] .+ spac_750.canopy.sun_geometry.auxil.e_difꜜ[:,irt] .+ spac_750.canopy.sun_geometry.auxil.e_difꜛ[:,irt];
        rad_cyto = spac_750_cyto.canopy.sun_geometry.auxil.e_dirꜜ[:,irt] .+ spac_750_cyto.canopy.sun_geometry.auxil.e_difꜜ[:,irt] .+ spac_750_cyto.canopy.sun_geometry.auxil.e_difꜛ[:,irt];
        leaf_750 = spac_750.plant.leaves[ilf];
        leaf_700 = spac_700.plant.leaves[ilf];
        leaf_750_cyto = spac_750_cyto.plant.leaves[ilf];
        leaf_700_cyto = spac_700_cyto.plant.leaves[ilf];
        f_sunlit = spac_750.canopy.sun_geometry.s_aux.p_sunlit[irt];
        push!(df, [sum(rad[mask_uvio]) / 1000,
                   sum(rad[mask_rpar]) / 1000,
                   sum(rad[mask_farr]) / 1000,
                   sum(rad[mask_grin]) / 1000,
                   mean(leaf_750.flux.auxil.a_g_sunlit) * f_sunlit + leaf_750.flux.auxil.a_g_shaded * (1 - f_sunlit),
                   mean(leaf_700.flux.auxil.a_g_sunlit) * f_sunlit + leaf_700.flux.auxil.a_g_shaded * (1 - f_sunlit),
                   mean(leaf_750.flux.state.g_H₂O_s_sunlit) * f_sunlit + leaf_750.flux.state.g_H₂O_s_shaded * (1 - f_sunlit),
                   mean(leaf_700.flux.state.g_H₂O_s_sunlit) * f_sunlit + leaf_700.flux.state.g_H₂O_s_shaded * (1 - f_sunlit),
                   sum(rad_cyto[mask_uvio]) / 1000,
                   sum(rad_cyto[mask_rpar]) / 1000,
                   sum(rad_cyto[mask_farr]) / 1000,
                   sum(rad_cyto[mask_grin]) / 1000,
                   mean(leaf_750_cyto.flux.auxil.a_g_sunlit) * f_sunlit + leaf_750_cyto.flux.auxil.a_g_shaded * (1 - f_sunlit),
                   mean(leaf_700_cyto.flux.auxil.a_g_sunlit) * f_sunlit + leaf_700_cyto.flux.auxil.a_g_shaded * (1 - f_sunlit),
                   mean(leaf_750_cyto.flux.state.g_H₂O_s_sunlit) * f_sunlit + leaf_750_cyto.flux.state.g_H₂O_s_shaded * (1 - f_sunlit),
                   mean(leaf_700_cyto.flux.state.g_H₂O_s_sunlit) * f_sunlit + leaf_700_cyto.flux.state.g_H₂O_s_shaded * (1 - f_sunlit)]);
    end;

    # 2.3 save the data
    save_csv!("../output/5_fr_diff_profiles.csv", df);
end;
