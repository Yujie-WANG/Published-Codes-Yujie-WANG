using CanopyLayers: FourBandsFittingCurve, FourBandsFittingHybrid, FourBandsFittingPoint, SoilOpticals, TwoBandsFittingCurve, TwoBandsFittingHybrid, TwoBandsFittingPoint
using CanopyLayers: create_wave_length, fit_soil_mat!, soil_albedos
using PkgUtility: nanmean, nanstd, read_csv, rmse
using PlotPlants: create_canvas, plot_line_regress!, save_canvas!, set_xlabels!, set_xticks!, set_xticklabels!, set_xylabels!, set_xylims!, set_xyticks!, set_ylabels!, set_ylims!, use_serif_tex


COLORS   = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"];


"""
    plot_SA_methods!(;saving::Bool = false, use_latex::Bool = true)

Plot the soil albedo project figures, given
- `proj` [`SoilAlbedo2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_SA_methods!(FT = Float64; saving::Bool = false, use_latex::Bool = true)
    # use latex and serif
    if use_latex use_serif_tex(); end

    # create soil optical structures for fitting purpose
    _wls = create_wave_length(FT);
    _soil_2p = SoilOpticals{FT}(_wls);
    _soil_2c = SoilOpticals{FT}(_wls);
    _soil_2h = SoilOpticals{FT}(_wls);
    _soil_4p = SoilOpticals{FT}(_wls);
    _soil_4c = SoilOpticals{FT}(_wls);
    _soil_4h = SoilOpticals{FT}(_wls);

    #
    #
    # raw bands from GSV model
    #
    #
    # plot the raw bands and save them as a figure
    _fig,(_ax1,) = create_canvas("bands");
    _ax1.plot(400:10:2500, _soil_2p.SW_mat_raw_4[:,1], label="Band-1");
    _ax1.plot(400:10:2500, _soil_2p.SW_mat_raw_4[:,2], label="Band-2");
    _ax1.plot(400:10:2500, _soil_2p.SW_mat_raw_4[:,3], label="Band-3");
    _ax1.plot(400:10:2500, _soil_2p.SW_mat_raw_4[:,4], label="Band-4");
    _ax1.legend(loc="lower right", ncol=2);

    set_ylims!([_ax1], [-1.15,1.05]);
    set_xylabels!([_ax1], "Wavelength (nm)", "PCA Bands");
    save_canvas!(_fig, "1-bands.pdf", saving);
    if saving close(_fig) end;

    #
    #
    # interpolation is working for all the curves
    #
    #
    # read the curves: 23 curves in total
    _df = read_csv("$(@__DIR__)/collections.csv");
    _sdfs = [];
    for _i in 1:23
        _sdf = _df[_df.ID .== _i,:];
        push!(_sdfs, _sdf);
    end;

    # interpolate the curves
    _refs = ones(FT, (length(_wls.WL),23)) .* FT(NaN);
    for _i in eachindex(_sdfs)
        for _j in eachindex(_wls.WL)
            _k = findfirst( (_sdfs[_i].WL .>= _wls.WL[_j]) );
            if !isnothing(_k) && _k > 1
                _x1 = _sdfs[_i].WL[_k];
                _x0 = _sdfs[_i].WL[_k-1];
                _y1 = _sdfs[_i].Albedo[_k];
                _y0 = _sdfs[_i].Albedo[_k-1];
                _x  = _wls.WL[_j];
                _al = (_x1 - _x) / (_x1 - _x0) * _y0 + (_x - _x0) / (_x1 - _x0) * _y1;
                _refs[_j,_i] = _al;
            end;
        end;
    end;

    _fig,(_ax1,) = create_canvas("curves"; figsize=(6,6));
    for _i in eachindex(_sdfs)
        _ax1.plot(_sdfs[_i].WL, _sdfs[_i].Albedo, "-"; alpha=0.5);
        _ax1.plot(_wls.WL, _refs[:,_i], "k:");
    end;

    set_xylabels!([_ax1], "Wavelength (nm)", "Albedo (-)");
    save_canvas!(_fig, "2-interpolation.pdf", saving);
    if saving close(_fig) end;

    #
    #
    # how well the methods work
    #
    #
    _rmse_2p = FT[];
    _rmse_2c = FT[];
    _rmse_2h = FT[];
    _rmse_4p = FT[];
    _rmse_4c = FT[];
    _rmse_4h = FT[];

    # plot the fitting per data using raw curves
    _fig,_axs = create_canvas("fittings"; nrow=4, ncol=6, figsize=(15,9));
    for _i in eachindex(_sdfs)
        _sdf = _sdfs[_i];
        _soil_2p.color = _sdf.Color[1];
        _soil_2c.color = _sdf.Color[1];
        _soil_2h.color = _sdf.Color[1];
        _soil_4p.color = _sdf.Color[1];
        _soil_4c.color = _sdf.Color[1];
        _soil_4h.color = _sdf.Color[1];
        _ref_nir = nanmean(_refs[_wls.iNIR,_i]);
        _ref_par = nanmean(_refs[_wls.iPAR,_i]);
        fit_soil_mat!(_soil_2p, _wls, _ref_par, _ref_nir, TwoBandsFittingPoint());
        fit_soil_mat!(_soil_2c, _wls, _ref_par, _ref_nir, TwoBandsFittingCurve());
        fit_soil_mat!(_soil_2h, _wls, _ref_par, _ref_nir, TwoBandsFittingHybrid());
        fit_soil_mat!(_soil_4p, _wls, _ref_par, _ref_nir, FourBandsFittingPoint());
        fit_soil_mat!(_soil_4c, _wls, _ref_par, _ref_nir, FourBandsFittingCurve());
        fit_soil_mat!(_soil_4h, _wls, _ref_par, _ref_nir, FourBandsFittingHybrid());
        _axs[_i].plot(_wls.WL, _refs[:,_i], "k:");
        _axs[_i].plot(_wls.WL, _soil_2p.ρ_SW, "-", alpha=0.7, label="2P");
        _axs[_i].plot(_wls.WL, _soil_2c.ρ_SW, "-", alpha=0.7, label="2C");
        _axs[_i].plot(_wls.WL, _soil_2h.ρ_SW, "-", alpha=0.7, label="2H");
        _axs[_i].plot(_wls.WL, _soil_4p.ρ_SW, "-", alpha=0.3, label="4P");
        _axs[_i].plot(_wls.WL, _soil_4c.ρ_SW, "-", alpha=0.7, label="4C");
        _axs[_i].plot(_wls.WL, _soil_4h.ρ_SW, "-", alpha=0.7, label="4H");
        push!(_rmse_2p, rmse(_refs[:,_i], _soil_2p.ρ_SW));
        push!(_rmse_2c, rmse(_refs[:,_i], _soil_2c.ρ_SW));
        push!(_rmse_2h, rmse(_refs[:,_i], _soil_2h.ρ_SW));
        push!(_rmse_4p, rmse(_refs[:,_i], _soil_4p.ρ_SW));
        push!(_rmse_4c, rmse(_refs[:,_i], _soil_4c.ρ_SW));
        push!(_rmse_4h, rmse(_refs[:,_i], _soil_4h.ρ_SW));
    end;
    _axs[end].plot(0, 0, "k:", label="obs");
    _axs[end].plot(0, 0, "-", alpha=0.7, label="2P");
    _axs[end].plot(0, 0, "-", alpha=0.7, label="2C");
    _axs[end].plot(0, 0, "-", alpha=0.7, label="2H");
    _axs[end].plot(0, 0, "-", alpha=0.3, label="4P");
    _axs[end].plot(0, 0, "-", alpha=0.7, label="4C");
    _axs[end].plot(0, 0, "-", alpha=0.7, label="4H");
    _axs[end].spines["bottom"].set_visible(false);
    _axs[end].spines["left"  ].set_visible(false);
    _axs[end].spines["right" ].set_visible(false);
    _axs[end].spines["top"   ].set_visible(false);
    set_xyticks!([_axs[end]], [[]], [[]]);
    _axs[end].legend(loc="center");

    set_xlabels!(_axs[19:end-1], "Wavelength (nm)");
    set_ylabels!(_axs[[1,7,13,19]], "Albedo (-)");
    set_ylims!(_axs, [0,0.5]);
    save_canvas!(_fig, "3-fitting.pdf", saving);
    if saving close(_fig) end;

    #
    #
    # the mean and std of RMSE for 6 methods
    #
    #
    _rmses = [_rmse_2p, _rmse_2c, _rmse_2h, _rmse_4p, _rmse_4c, _rmse_4h];
    _fig,(_ax1,) = create_canvas("fitting errors");
    _ax1.bar(1:6, [nanmean(_dat) for _dat in _rmses], color=COLORS[1:6], yerr=[nanstd(_dat) for _dat in _rmses],  alpha=0.7);

    set_xticks!([_ax1], [1,2,3,4,5,6]);
    set_xticklabels!([_ax1], ["2P", "2C", "2H", "4P", "4C", "4H"]);
    set_ylabels!([_ax1], "RMSE of Albedo (-)");
    set_ylims!([_ax1], [0,0.06]);
    save_canvas!(_fig, "4-fitting-error.pdf", saving);
    if saving close(_fig) end;

    #
    #
    # how well the CLM table works
    #
    #
    _ρ_par_obs = FT[];
    _ρ_nir_obs = FT[];
    _ρ_par_clm = FT[];
    _ρ_nir_clm = FT[];
    _swc_label = FT[];
    for _i in eachindex(_sdfs)
        _sdf = _sdfs[_i];
        push!(_ρ_nir_obs, nanmean(_refs[_wls.iNIR,_i]));
        push!(_ρ_par_obs, nanmean(_refs[_wls.iPAR,_i]));
        _clm_ρs = soil_albedos(_sdf.Color[1], _sdf.VWC[1], true);
        push!(_ρ_nir_clm, _clm_ρs[2]);
        push!(_ρ_par_clm, _clm_ρs[1]);
        push!(_swc_label, _sdfs[_i].VWC[1]);
    end;

    _colors = [1,1,1,2,2,2,2,2,3,3,3,3,3,3,4,4,4,5,5,5,5,5,5];

    _fig,(_ax1,) = create_canvas("broadband errors");
    for _i in eachindex(_ρ_nir_clm)
        _ax1.plot(_ρ_par_obs[_i], _ρ_par_clm[_i], "o", color=COLORS[_colors[_i]], alpha=0.2+0.8*_swc_label[_i]);
        _ax1.plot(_ρ_nir_obs[_i], _ρ_nir_clm[_i], "s", color=COLORS[_colors[_i]], alpha=0.2+0.8*_swc_label[_i]);
    end;
    _ax1.plot(-0.1,-0.1, "ko", label="PAR");
    _ax1.plot(-0.1,-0.1, "ks", label="NIR");
    _ax1.plot([0,0.5], [0,0.5], "k-");
    plot_line_regress!(_ax1, [_ρ_par_obs; _ρ_nir_obs], [_ρ_par_clm; _ρ_nir_clm], color="gray", interval=true);
    _ax1.legend(loc="lower right");

    set_xylims!([_ax1], [0,0.5], [0,0.5]);
    set_xylabels!([_ax1], "obs Albedo (-)", "CliMA Albedo (-)");
    save_canvas!(_fig, "5-clm.pdf", saving);
    if saving close(_fig) end;

    #
    #
    # how well the methods work directly using CLM
    #
    #
    _rmse_2p = FT[];
    _rmse_2c = FT[];
    _rmse_2h = FT[];
    _rmse_4p = FT[];
    _rmse_4c = FT[];
    _rmse_4h = FT[];

    # plot the fitting per data using raw curves
    _fig,_axs = create_canvas("clm"; nrow=3, ncol=7, figsize=(18,9));
    _slection = [collect(1:14); collect(18:23)];
    @show _slection;
    for _i in eachindex(_slection)
        _sdf = _sdfs[_slection[_i]];
        _soil_2p.color = _sdf.Color[1];
        _soil_2c.color = _sdf.Color[1];
        _soil_2h.color = _sdf.Color[1];
        _soil_4p.color = _sdf.Color[1];
        _soil_4c.color = _sdf.Color[1];
        _soil_4h.color = _sdf.Color[1];
        fit_soil_mat!(_soil_2p, _wls, _sdf.VWC[1], TwoBandsFittingPoint());
        fit_soil_mat!(_soil_2c, _wls, _sdf.VWC[1], TwoBandsFittingCurve());
        fit_soil_mat!(_soil_2h, _wls, _sdf.VWC[1], TwoBandsFittingHybrid());
        fit_soil_mat!(_soil_4p, _wls, _sdf.VWC[1], FourBandsFittingPoint());
        fit_soil_mat!(_soil_4c, _wls, _sdf.VWC[1], FourBandsFittingCurve());
        fit_soil_mat!(_soil_4h, _wls, _sdf.VWC[1], FourBandsFittingHybrid());
        _axs[_i].plot(_wls.WL, _refs[:,_slection[_i]], "k:");
        _axs[_i].plot(_wls.WL, _soil_2p.ρ_SW, "-", alpha=0.7, label="2P");
        _axs[_i].plot(_wls.WL, _soil_2c.ρ_SW, "-", alpha=0.7, label="2C");
        _axs[_i].plot(_wls.WL, _soil_2h.ρ_SW, "-", alpha=0.7, label="2H");
        _axs[_i].plot(_wls.WL, _soil_4p.ρ_SW, "-", alpha=0.3, label="4P");
        _axs[_i].plot(_wls.WL, _soil_4c.ρ_SW, "-", alpha=0.7, label="4C");
        _axs[_i].plot(_wls.WL, _soil_4h.ρ_SW, "-", alpha=0.7, label="4H");
        push!(_rmse_2p, rmse(_refs[:,_slection[_i]], _soil_2p.ρ_SW));
        push!(_rmse_2c, rmse(_refs[:,_slection[_i]], _soil_2c.ρ_SW));
        push!(_rmse_2h, rmse(_refs[:,_slection[_i]], _soil_2h.ρ_SW));
        push!(_rmse_4p, rmse(_refs[:,_slection[_i]], _soil_4p.ρ_SW));
        push!(_rmse_4c, rmse(_refs[:,_slection[_i]], _soil_4c.ρ_SW));
        push!(_rmse_4h, rmse(_refs[:,_slection[_i]], _soil_4h.ρ_SW));
    end;
    _axs[end].plot(0, 0, "k:", label="obs");
    _axs[end].plot(0, 0, "-", alpha=0.7, label="2P");
    _axs[end].plot(0, 0, "-", alpha=0.7, label="2C");
    _axs[end].plot(0, 0, "-", alpha=0.7, label="2H");
    _axs[end].plot(0, 0, "-", alpha=0.3, label="4P");
    _axs[end].plot(0, 0, "-", alpha=0.7, label="4C");
    _axs[end].plot(0, 0, "-", alpha=0.7, label="4H");
    _axs[end].spines["bottom"].set_visible(false);
    _axs[end].spines["left"  ].set_visible(false);
    _axs[end].spines["right" ].set_visible(false);
    _axs[end].spines["top"   ].set_visible(false);
    set_xyticks!([_axs[end]], [[]], [[]]);
    _axs[end].legend(loc="center");

    set_xlabels!(_axs[15:end-1], "Wavelength (nm)");
    set_ylabels!(_axs[[1,8,15]], "Albedo (-)");
    set_ylims!(_axs, [0,0.5]);
    save_canvas!(_fig, "6-clm.pdf", saving);
    if saving close(_fig) end;

    #
    #
    # the mean and std of RMSE for 6 methods
    #
    #
    _rmses = [_rmse_2p, _rmse_2c, _rmse_2h, _rmse_4p, _rmse_4c, _rmse_4h];
    _fig,(_ax1,) = create_canvas("clm errors");
    _ax1.bar(1:6, [nanmean(_dat) for _dat in _rmses], color=COLORS[1:6], yerr=[nanstd(_dat) for _dat in _rmses],  alpha=0.7);

    set_xticks!([_ax1], [1,2,3,4,5,6]);
    set_xticklabels!([_ax1], ["2P", "2C", "2H", "4P", "4C", "4H"]);
    set_ylabels!([_ax1], "RMSE of Albedo (-)");
    set_ylims!([_ax1], [0,0.1]);
    save_canvas!(_fig, "7-clm-error.pdf", saving);
    if saving close(_fig) end;

    return nothing
end




function fit_SA_methods2!(FT = Float64; saving::Bool = false, use_latex::Bool = true)
    # use latex and serif
    if use_latex use_serif_tex(); end

    # create soil optical structures for fitting purpose
    _wls = create_wave_length(FT);

    # read the curves: 23 curves in total
    _df = read_csv("$(@__DIR__)/collections.csv");
    _sdfs = [];
    for _i in 1:23
        _sdf = _df[_df.ID .== _i,:];
        push!(_sdfs, _sdf);
    end;

    # interpolate the curves
    _refs = ones(FT, (length(_wls.WL),23)) .* FT(NaN);
    for _i in eachindex(_sdfs)
        for _j in eachindex(_wls.WL)
            _k = findfirst( (_sdfs[_i].WL .>= _wls.WL[_j]) );
            if !isnothing(_k) && _k > 1
                _x1 = _sdfs[_i].WL[_k];
                _x0 = _sdfs[_i].WL[_k-1];
                _y1 = _sdfs[_i].Albedo[_k];
                _y0 = _sdfs[_i].Albedo[_k-1];
                _x  = _wls.WL[_j];
                _al = (_x1 - _x) / (_x1 - _x0) * _y0 + (_x - _x0) / (_x1 - _x0) * _y1;
                _refs[_j,_i] = _al;
            end;
        end;
    end;

    # how well the fitted table works
    _ρ_par_obs = FT[];
    _ρ_nir_obs = FT[];
    _ρ_par_clm = FT[];
    _ρ_nir_clm = FT[];
    _swc_label = FT[];
    for _i in eachindex(_sdfs)
        _sdf = _sdfs[_i];
        push!(_ρ_nir_obs, nanmean(_refs[_wls.iNIR,_i]));
        push!(_ρ_par_obs, nanmean(_refs[_wls.iPAR,_i]));
        _clm_ρs = soil_albedos(_sdf.Color[1], _sdf.VWC[1]);
        push!(_ρ_nir_clm, _clm_ρs[2]);
        push!(_ρ_par_clm, _clm_ρs[1]);
        push!(_swc_label, _sdfs[_i].VWC[1]);
    end;

    _colors = [1,1,1,2,2,2,2,2,3,3,3,3,3,3,4,4,4,5,5,5,5,5,5];

    _fig,(_ax1,) = create_canvas("x broadband errors");
    for _i in eachindex(_ρ_nir_clm)
        _ax1.plot(_ρ_par_obs[_i], _ρ_par_clm[_i], "o", color=COLORS[_colors[_i]], alpha=0.2+0.8*_swc_label[_i]);
        _ax1.plot(_ρ_nir_obs[_i], _ρ_nir_clm[_i], "s", color=COLORS[_colors[_i]], alpha=0.2+0.8*_swc_label[_i]);
    end;
    _ax1.plot(-0.1,-0.1, "ko", label="PAR");
    _ax1.plot(-0.1,-0.1, "ks", label="NIR");
    _ax1.plot([0,0.5], [0,0.5], "k-");
    plot_line_regress!(_ax1, [_ρ_par_obs; _ρ_nir_obs], [_ρ_par_clm; _ρ_nir_clm], color="gray", interval=true);
    _ax1.legend(loc="lower right");

    set_xylims!([_ax1], [0,0.5], [0,0.5]);
    set_xylabels!([_ax1], "obs Albedo (-)", "CLM Albedo (-)");
    save_canvas!(_fig, "x-clm.pdf", saving);
    if saving close(_fig) end;

    return nothing
end
