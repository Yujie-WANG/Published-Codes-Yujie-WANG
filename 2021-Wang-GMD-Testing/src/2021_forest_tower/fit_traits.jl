###############################################################################
#
# Threading functions to fit Vcmax, Rbase, Kmax
#
###############################################################################
"""
    fit_vrk!(proj::ForestTower2021{FT},
             node::SPACMono{FT},
             df::DataFrame,
             year::Int,
             debugging::Bool,
             site::String,
             label::String
    ) where {FT<:AbstractFloat}

Fit Vcmax, Rbase, and Kmax, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `node` `SPACMono` type struct
- `df` DataFrame that contains data and output
- `year` Which year
- `debugging` If true, dispay the output per vrk combination
- `site` Site of flux tower, must be `NiwotRidge` or `Ozark`
- `label` Label to show up on generated CSV files to distinct stomatal models

Note here the fitted `Kmax` is twice the whole plant maximum.
"""
function fit_vrk!(
            proj::ForestTower2021{FT},
            node::SPACMono{FT},
            df::DataFrame,
            year::Int,
            debugging::Bool,
            site::String,
            label::String
) where {FT<:AbstractFloat}
    # use 1.4 for angiosperm and 1.7 for gymnosperm
    _q10 = (site=="NiwotRidge" ? 1.7 : 1.4);

    # minimize the error to infer leaf and soil respiration
    @inline f(x) = (
        # Kmax fitted here is the root Kmax
        update_VJRWW!(node, x[1]);
        update_Kmax!(node, x[3] / 2);
        rbase = Q10TD{FT}(x[2], 298.15, _q10);

        # ignore NaN values from the diff series
        simulation!(proj, deepcopy(node), df, site, rbase);
        _diff_c = abs.( (df.ObsC .- df.ModC) ./ nanstd(df.ObsC) );
        _diff_e = abs.( (df.ObsE .- df.ModE) ./ nanstd(df.ObsE) );
        _err    = sum(_diff_c[.~isnan.(_diff_c)]) +
                  sum(_diff_e[.~isnan.(_diff_e)]);

        if debugging
            @show x,_err;
        end;
        # uncomment this to show the fitting process
        #figure(1); clf();
        #subplot(2,1,1); plot(df.ObsC); plot(df.ModC);
        #subplot(2,1,2); plot(df.ObsE); plot(df.ModE);

        return -_err
    );

    # fit the data
    ms  = ReduceStepMethodND{FT}(
                 x_mins = FT[15, 1, 0.02],
                 x_maxs = FT[100, 100, 8],
                 x_inis = FT[20, 2, 0.2],
                 Δ_inis = FT[10, 10, 10]);
    st  = SolutionToleranceND{FT}(FT[0.101, 0.101, 0.101], 50);
    sol = find_peak(f, ms, st);

    return [site label year sol[1] sol[2] sol[3] NaN], df.ModC, df.ModE
end








###############################################################################
#
# Threading functions to fit Vcmax, Rbase, Kmax, and g1
#
###############################################################################
"""
    fit_vrkg!(proj::ForestTower2021{FT},
              node::SPACMono{FT},
              df::DataFrame,
              year::Int,
              debugging::Bool,
              site::String,
              label::String
    ) where {FT<:AbstractFloat}

Fit Vcmax, Rbase, Kmax, and g1 (for Ball-Berry and Medlyn models only), given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `node` `SPACMono` type struct
- `df` DataFrame that contains data and output
- `year` Which year
- `debugging` If true, dispay the output per vrk combination
- `site` Site of flux tower, must be `NiwotRidge` or `Ozark`
- `label` Label to show up on generated CSV files to distinct stomatal models

Note here the fitted `Kmax` is twice the whole plant maximum.
"""
function fit_vrkg!(
            proj::ForestTower2021{FT},
            node::SPACMono{FT},
            df::DataFrame,
            year::Int,
            debugging::Bool,
            site::String,
            label::String
) where {FT<:AbstractFloat}
    # use 1.4 for angiosperm and 1.7 for gymnosperm
    _q10 = (site=="NiwotRidge" ? 1.7 : 1.4);

    # minimize the error to infer leaf and soil respiration
    @inline f(x) = (
        # Kmax fitted here is root Kmax
        update_VJRWW!(node, x[1]);
        update_Kmax!(node, x[3] / 2);
        rbase = Q10TD{FT}(x[2], 298.15, _q10);

        # g1 for empirical stomatal model
        node.stomata_model.g1 = x[4];

        # ignore NaN values from the diff series
        simulation!(proj, deepcopy(node), df, site, rbase);
        _diff_c = abs.( (df.ObsC .- df.ModC) ./ nanstd(df.ObsC) );
        _diff_e = abs.( (df.ObsE .- df.ModE) ./ nanstd(df.ObsE) );
        _err    = sum(_diff_c[.~isnan.(_diff_c)]) +
                  sum(_diff_e[.~isnan.(_diff_e)]);

        if debugging
            @show x,_err;
        end;
        # uncomment this to show the fitting process
        #figure(1); clf();
        #subplot(2,1,1); plot(df.ObsC); plot(df.ModC);
        #subplot(2,1,2); plot(df.ObsE); plot(df.ModE);

        return -_err
    );

    # fit the data for different models
    if typeof(node.stomata_model) <: ESMMedlyn
        g_ini = 200;
        g_max = 1000;
        g_del = 100;
        g_lim = 1.01;
    else
        g_ini = 20;
        g_max = 100;
        g_del = 10;
        g_lim = 0.101;
    end
    ms  = ReduceStepMethodND{FT}(
                 x_mins = FT[15, 1, 0.02, 1],
                 x_maxs = FT[100, 100, 8, g_max],
                 x_inis = FT[20, 2, 0.2, g_ini],
                 Δ_inis = FT[10, 10, 10, g_del]);
    st  = SolutionToleranceND{FT}(FT[0.101, 0.101, 0.101, g_lim], 30);
    sol = find_peak(f, ms, st);

    return [site label year sol[1] sol[2] sol[3] sol[4]], df.ModC, df.ModE
end








###############################################################################
#
# Threading functions to fit Vcmax, Rbase, Kmax, and g1 for all setups
#
###############################################################################
"""
Fit the data for flux towers.
"""
function fit_tower! end




"""
This method is a wrapper for [`fit_vrk!`](@ref) and [`fit_vrkg!`](@ref). If the
    7th item of the input vector is `bbmg` or `medg`, [`fit_vrkg!`](@ref) is
    used; otherwise, [`fit_vrk!`](@ref) is used

    fit_tower!(param::Array)

A wrapper to parse parameters to threading functions, given
- `params` A vector of parameters
"""
fit_tower!(param::Array) =
(
    #proj, deepcopy(node), data, year, debugging, site, label
    @assert param[6] in ["NiwotRidge", "Ozark"];
    @assert param[7] in ["bbm", "bbmg", "med", "medg", "osm"];

    if param[7] in ["bbmg", "medg"]
        return fit_vrkg!(param...)
    else
        return fit_vrk!(param...)
    end;
)




"""
This method runs all the simulations in parallel and save the simulation
    results as two files: one for fitting results, and one for simulation time
    series. Each model setup is labeled by model name and label.

    fit_tower!(proj::ForestTower2021{FT},
               nTH::Int = 130,
               debugging::Bool = true
    ) where {FT<:AbstractFloat}

Fit all the flux data at the same time using multi-threading calculation, given
- `proj` A [`ForestTower2021`](@ref) type project indentifier
- `nTH` Number of threadings to run in parallel. Default is `130`.
- `debugging` If true, run the model in debugging mode. Default is `true`.
"""
fit_tower!(proj::ForestTower2021{FT},
           nTH::Int = 130,
           debugging::Bool = true
) where {FT<:AbstractFloat} =
(
    @info tinfo("Loading $(nTH) threadings...");
    dynamic_workers!(proj, nTH);

    # create parameter list
    _params = query_data(proj, debugging);
    if debugging
        _fittings = pmap(fit_tower!, _params);
    else
        _fittings = @showprogress pmap(fit_tower!, _params);
    end;

    # parse results
    save_csv!("fitted_all.csv", [_r[1] for _r in _fittings],
              ["Site", "Model", "Year", "Vcmax", "Rbase", "K", "g1"]);
    for _i in eachindex(_params)
        _params[_i][3].ModC = _fittings[_i][2];
        _params[_i][3].ModE = _fittings[_i][3];
    end;
    _all_df = vcat( [_para[3] for _para in _params]... );
    save_csv!("simulated_all.csv", _all_df);

    # send out a notification email
    send_email!("Project Status", "fluo@gps.caltech.edu", "jesiner@gmail.com",
                "Function fit_tower! simulations done!");

    return nothing
)








###############################################################################
#
# Reprocess CSV files from all results to individual files
#
###############################################################################
"""
    divide_results!(proj::ForestTower2021{FT}) where {FT<:AbstractFloat}

Divide the simulated results to multiple files so as to be used for plotting
    purpose, given
- `proj` [`ForestTower2021`](@ref) type project control

Note that the simulations are done on the server, so remember to copy the newly
    simulated results to local PC before running this function.
"""
function divide_results!(proj::ForestTower2021{FT}) where {FT<:AbstractFloat}
    @info tinfo("Divide the simulation results to parts...");

    # read files
    _folder = "/home/wyujie/RAID/Data/FLUXNET2015/simulation/";
    _data_fit = read_csv("$(_folder)fitted_all.csv");
    _data_all = read_csv("$(_folder)simulated_all.csv");

    # function to save subsets
    @inline save_subset!(site, model, lab1, lab2) = (
        _mask_fit = (_data_fit.Site .== site) .* (_data_fit.Model .== model);
        _mask_all = (_data_all.Site .== site) .* (_data_all.Label .== model);
        _set_fit = _data_fit[_mask_fit, :];
        _set_all = _data_all[_mask_all, :];

        save_csv!(_set_fit, "$(_folder)fitted_$(lab1)_$(site).csv");
        save_csv!(_set_all, "$(_folder)simulated_results_$(lab2)_$(site).csv");

        return nothing;
    );

    # save the data per model set up
    save_subset!("NiwotRidge", "bbm", "vrk_bbm", "bbm");
    save_subset!("NiwotRidge", "med", "vrk_med", "med");
    save_subset!("NiwotRidge", "osm", "vrk_osm", "osm");
    save_subset!("NiwotRidge", "bbmg", "vrkg_bbm", "g_bbm");
    save_subset!("NiwotRidge", "medg", "vrkg_med", "g_med");
    save_subset!("Ozark", "bbm", "vrk_bbm", "bbm");
    save_subset!("Ozark", "med", "vrk_med", "med");
    save_subset!("Ozark", "osm", "vrk_osm", "osm");
    save_subset!("Ozark", "bbmg", "vrkg_bbm", "g_bbm");
    save_subset!("Ozark", "medg", "vrkg_med", "g_med");

    return nothing
end
