#
# this script is meant to fit the SPAC to the data collected from Zhen et al. (2020) PCE
#
using Distributed: @everywhere, addprocs, pmap

using Emerald.EmeraldUtility.Threading: dynamic_workers!

output_data = "../output/6_fitting.csv";
if !isfile(output_data)
    dynamic_workers!(32);
end;
@everywhere FT = Float64;
@everywhere using DataFrames: DataFrame
@everywhere using Statistics: mean
@everywhere using Emerald.EmeraldIO.Text: read_csv, save_csv!
@everywhere using Emerald.EmeraldLand.Namespace: BulkSPAC, C3CytoState, C3CytoTrait, LAND_2021_1NM, SPACConfiguration
@everywhere using Emerald.EmeraldLand.SPAC: GPP, PPAR, initialize_spac!, prescribe_air!, prescribe_traits!, read_spectrum, soil_plant_air_continuum!
@everywhere using Emerald.EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak
@everywhere using Emerald.EmeraldMath.Stats: rmse
@everywhere using Emerald.EmeraldPhysics.Optics: photon


# 1. read the white light and far red light spectra and interpolate them to the model's wavelength per 1 nm
@everywhere spectra_file = "../output/6_spectra.csv";
@everywhere CONFIG = SPACConfiguration(FT, dataset = LAND_2021_1NM);
@everywhere CONFIG.ENABLE_REF = false;
@everywhere CONFIG.ENABLE_SIF = false;
@everywhere spectrum_par = similar(CONFIG.SPECTRA.Λ);
@everywhere spectrum_far = similar(CONFIG.SPECTRA.Λ);
local_par = similar(CONFIG.SPECTRA.Λ);
local_far = similar(CONFIG.SPECTRA.Λ);
if !isfile(spectra_file)
    white_light = read_csv("../../../data/white-light.csv");
    farred_light = read_csv("../../../data/farred-light.csv");
    max_par,min_par = maximum(white_light.WL),minimum(white_light.WL);
    max_far,min_far = maximum(farred_light.WL),minimum(farred_light.WL);
    for i in eachindex(CONFIG.SPECTRA.Λ)
        if min_par <= CONFIG.SPECTRA.Λ[i] <= max_par
            local_par[i] = read_spectrum(white_light.WL, white_light.RAD, CONFIG.SPECTRA.Λ[i]);
        else
            local_par[i] = 0;
        end;
        if min_far <= CONFIG.SPECTRA.Λ[i] <= max_far
            local_far[i] = read_spectrum(farred_light.WL, farred_light.RAD, CONFIG.SPECTRA.Λ[i]);
        else
            local_far[i] = 0;
        end;
    end;
    save_csv!(spectra_file, DataFrame(WL = CONFIG.SPECTRA.Λ, PAR = local_par, FR = local_far));
end;
@everywhere spectra = read_csv(spectra_file);
@everywhere spectrum_par = spectra.PAR;
@everywhere spectrum_far = spectra.FR;
@everywhere photon_par = sum(photon.(CONFIG.SPECTRA.Λ, spectrum_par)) / 1000 * 1e6;
@everywhere photon_far = sum(photon.(CONFIG.SPECTRA.Λ, spectrum_far)) / 1000 * 1e6;

# 2. define function to compute the RMSE between the data and the SPAC model
@everywhere function GPPS(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, treedata::DataFrame) where {FT}
    # per PPFD, rescale the input radiance to the diffuse radiation in the meteo
    gpps = [];
    for i in eachindex(treedata.TPAR)
        tmp_spac = deepcopy(spac);
        tmp_spac.meteo.rad_sw.e_dir .= eps(FT);
        tmp_spac.meteo.rad_sw.e_dif .= treedata.PAR[i] / photon_par * spectrum_par .+ treedata.FR[i] / photon_far * spectrum_far;
        sw_rad = sum(tmp_spac.meteo.rad_sw.e_dir + tmp_spac.meteo.rad_sw.e_dif) / 1000;
        tmp_spac.meteo.rad_lw += 400 - sw_rad;
        soil_plant_air_continuum!(config, tmp_spac, FT(3600));
        push!(gpps, GPP(tmp_spac));
    end;

    return gpps
end;

# 3. use the SPAC model to fit the data (fit Vcmax and Chl only)
if !isfile(output_data)
    input_data = "../../../data/par-fr-gpp.csv";
    indata = read_csv(input_data);

    # iterate through the data per TREE
    all_trees = [];
    for i in unique(indata.TREE)
        tree_data = indata[indata.TREE .== i, :];

        # fit the SPAC model to the data
        spac = BulkSPAC(CONFIG);
        for alayer in spac.airs
            prescribe_air!(alayer; f_CO₂ = 600);
        end;
        spac_cyto = deepcopy(spac);
        for l in spac_cyto.plant.leaves
            l.photosystem.trait = C3CytoTrait{FT}();
            l.photosystem.state = C3CytoState{FT}();
        end;
        prescribe_traits!(CONFIG, spac; lai = tree_data.LAI[1], sai = 0, vcmax = 80, vcmax_expo = 0.3);
        for l in spac.plant.leaves
            l.photosystem.trait.j_max25 = l.photosystem.trait.v_cmax25 * 2;
        end;
        prescribe_traits!(CONFIG, spac_cyto; lai = tree_data.LAI[1], sai = 0, vcmax = 80, vcmax_expo = 0.3);
        initialize_spac!(CONFIG, spac);
        initialize_spac!(CONFIG, spac_cyto);
        spac.meteo.rad_sw.e_dir .= eps(FT);
        spac.meteo.rad_sw.e_dif .= eps(FT);
        spac_cyto.meteo.rad_sw.e_dir .= eps(FT);
        spac_cyto.meteo.rad_sw.e_dif .= eps(FT);
        soil_plant_air_continuum!(CONFIG, spac, FT(3600));
        soil_plant_air_continuum!(CONFIG, spac_cyto, FT(3600));

        # create an array of PPFD per white and far red light
        tree_data[!, "PAR"] .= NaN;
        tree_data[!, "FR"] .= NaN;
        ppfd_ini = tree_data.TPAR[findfirst(tree_data.LIGHT .== "PAR")];
        for j in eachindex(tree_data.TPAR)
            if tree_data.LIGHT[j] == "FR"
                tree_data.PAR[j] = 0;
                tree_data.FR[j] = tree_data.TPAR[j];
            elseif tree_data.LIGHT[j] == "PAR"
                tree_data.PAR[j] = tree_data.TPAR[j];
                tree_data.FR[j] = 0;
            else
                tree_data.PAR[j] = ppfd_ini;
                tree_data.FR[j] = tree_data.TPAR[j] - ppfd_ini;
            end;
        end;
        push!(all_trees, [spac, spac_cyto, tree_data]);
    end;

    # define local function from target_function
    @everywhere fit_spac(param) = (
        spac = param[1];
        spac_cyto = param[2];
        treedata = param[3];
        sol_met = ReduceStepMethodND{FT}(x_mins = [1,1], x_maxs = [300, 200], x_inis = [80, 80], Δ_inis = [10, 10]);
        sol_tol = SolutionToleranceND{FT}([1, 1], 50);
        sol_met_cyto = deepcopy(sol_met);
        sol_tol = deepcopy(sol_tol);
        sol_fun(x) = (
            prescribe_traits!(CONFIG, spac; cab = x[2], car = x[2]/7, vcmax = x[1], vcmax_expo = 0.3);
            for l in spac.plant.leaves
                l.photosystem.trait.j_max25 = l.photosystem.trait.v_cmax25 * 2;
            end;
            initialize_spac!(CONFIG, spac);
            _rmse = -rmse(treedata.GPP, GPPS(CONFIG, spac, treedata));
            @show  x, _rmse;
            return _rmse
        );
        sol_fun_cyto(x) = (
            prescribe_traits!(CONFIG, spac_cyto; cab = x[2], car = x[2]/7, vcmax = x[1], vcmax_expo = 0.3);
            initialize_spac!(CONFIG, spac_cyto);
            _rmse = -rmse(treedata.GPP, GPPS(CONFIG, spac_cyto, treedata));
            @show  x, _rmse;
            return _rmse
        );
        sol = find_peak(sol_fun, sol_met, sol_tol);
        sol_cyto = find_peak(sol_fun_cyto, sol_met_cyto, sol_tol);

        prescribe_traits!(CONFIG, spac; cab = sol[2], car = sol[2]/7, vcmax = sol[1], vcmax_expo = 0.3);
        for l in spac.plant.leaves
            l.photosystem.trait.j_max25 = l.photosystem.trait.v_cmax25 * 2;
        end;
        prescribe_traits!(CONFIG, spac_cyto; cab = sol_cyto[2], car = sol_cyto[2]/7, vcmax = sol_cyto[1], vcmax_expo = 0.3);
        initialize_spac!(CONFIG, spac);
        initialize_spac!(CONFIG, spac_cyto);
        gpps = GPPS(CONFIG, spac, treedata);
        gpps_cyto = GPPS(CONFIG, spac_cyto, treedata);
        treedata[!, "CHL_FIT"] .= sol[2];
        treedata[!, "VCMAX_FIT"] .= sol[1];
        treedata[!, "GPP_FIT"] .= gpps;
        treedata[!, "CHL_FIT_CYTO"] .= sol_cyto[2];
        treedata[!, "VCMAX_FIT_CYTO"] .= sol_cyto[1];
        treedata[!, "GPP_FIT_CYTO"] .= gpps_cyto;

        return treedata
    );

    # fit the data in parallel and remove the workers when done
    results = pmap(fit_spac, all_trees);
    dynamic_workers!(0);

    # save the results
    save_csv!(output_data, vcat(results...));
end;
