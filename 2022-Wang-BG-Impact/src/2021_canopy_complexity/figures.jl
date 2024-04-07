###############################################################################
#
# Plot the main text figures for 2021_canopy_complexity project
#
###############################################################################
"""
    plot_CC2021!(
                proj::CanopyComplexity2021{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot CanopyComplexity2021 figures, given
- `proj` [`CanopyComplexity2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_CC2021!(
            proj::CanopyComplexity2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    plot_CC2021_NL!(proj; saving=saving, use_latex=use_latex);
    plot_CC2021_APAR!(proj; saving=saving, use_latex=use_latex);
    plot_CC2021_SA!(proj; saving=saving, use_latex=use_latex, gradients=false);
    plot_CC2021_SA!(proj; saving=saving, use_latex=use_latex, gradients=true);
    plot_CC2021_SA_err!(proj; saving=saving, use_latex=use_latex);
    diurnal_cycles!(proj; saving=saving, use_latex=use_latex);
    spectra_diff!(proj, "9-brdf-12", 12; saving=saving, use_latex=use_latex);
    spectra_diff!(proj, "10-brdf-16", 16; saving=saving, use_latex=use_latex);
    plot_CC2021_SIF_RAD!(proj; saving=saving, use_latex=use_latex);

    return nothing
end




"""
    plot_CC2021_NL!(
                proj::CanopyComplexity2021{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot non-linear photosynthesis responses to APAR or Vcmax, given
- `proj` [`CanopyComplexity2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_CC2021_NL!(
            proj::CanopyComplexity2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    # use latex and serif
    if use_latex use_serif_tex(); end

    #
    #
    # Non-linear response to APAR
    #
    #
    @info tinfo("Plot APAR and Vcmax response curves...");
    _canopy = CanopyLayer{FT}(n_leaf=1, g_max25=0.3);
    _envir  = AirLayer{FT}();
    _hs     = LeafHydraulics{FT}();
    _psm    = C3CLM(FT);
    _sm     = OSMWang{FT}();
    update_leaf_TP!(_psm, _canopy, _hs, _envir);
    temperature_effects!(_hs, _canopy.T);

    _apars = collect(FT, 10:10:500);
    _ans   = similar(_apars);
    _gsws  = similar(_apars);
    @showprogress for _i in eachindex(_apars)
        _canopy.APAR .= _apars[_i];
        _canopy.ps.APAR = _apars[_i];
        leaf_ETR!(_psm, _canopy.ps);
        for _n in 1:600
            gas_exchange!(_psm, _canopy, _envir, GswDrive());
            for _iLF in eachindex(_canopy.An)
                _a1 = _canopy.Ac[_iLF];
                _a2 = _canopy.Aj[_iLF];
                _a3 = _canopy.Ap[_iLF];
                _ag = lower_quadratic(0.99, -(_a1 + _a2), _a1 * _a2);
                _ag = lower_quadratic(0.99, -(_ag + _a3), _ag * _a3);
                _canopy.Ag[_iLF] = _ag;
                _canopy.An[_iLF] = _ag - _canopy.ps.Rd;
            end;
            update_gsw!(_canopy, _sm, _psm, _envir, FT(60); τ=FT(1e-6), smoothing=true);
            gsw_control!(_psm, _canopy, _envir);
        end;
        #gas_exchange!(_psm, _canopy, _hs, _envir, _sm);
        for _iLF in eachindex((_canopy).An)
            _a1 = (_canopy).Ac[_iLF];
            _a2 = (_canopy).Aj[_iLF];
            _a3 = (_canopy).Ap[_iLF];
            _ag = lower_quadratic(0.99, -(_a1 + _a2), _a1 * _a2);
            _ag = lower_quadratic(0.99, -(_ag + _a3), _ag * _a3);
            (_canopy).Ag[_iLF] = _ag;
            (_canopy).An[_iLF] = _ag - (_canopy).ps.Rd;
        end;
        _ans[_i]  = _canopy.An[1];
        _gsws[_i] = _canopy.g_sw[1];
    end

    # save figure
    _fig,_axs = create_canvas("CC-2"; figsize=(8,3.5), ncol=2);
    _ax1,_ax2 = _axs;
    _tx1 = _ax1.twinx();
    _ax1.plot(_apars, _ans, "g-", label=LS_Anet);
    _ax1.plot(0, 0, "c-", label="\$g_\\text{sw}\$");
    _tx1.plot(_apars, _gsws, "c-");
    _ax1.plot([50,50], [-1,22], "k:");
    _ax1.plot([350,350], [-1,22], "k:");

    _ax1.plot([50,200,350], [_ans[5],(_ans[5]+_ans[35])/2,_ans[35]], "g:");
    _ax1.plot([50,200,350], [_ans[5],(_ans[5]+_ans[35])/2,_ans[35]], "go");
    _ax1.plot(200, _ans[20], "go", mfc="none");

    _tx1.plot([50,200,350], [_gsws[5],(_gsws[5]+_gsws[35])/2,_gsws[35]], "c:");
    _tx1.plot([50,200,350], [_gsws[5],(_gsws[5]+_gsws[35])/2,_gsws[35]], "co");
    _tx1.plot(200, _gsws[20], "co", mfc="none");

    _ax1.annotate("", xy=(200,_ans[20]*0.98) , xytext=(200,(_ans[5]+_ans[35])/2*1.02)  , arrowprops=ARROW_G, fontsize=16, color="g");
    _tx1.annotate("", xy=(200,_gsws[20]*0.98), xytext=(200,(_gsws[5]+_gsws[35])/2*1.02), arrowprops=ARROW_C, fontsize=16, color="c");

    _ax1.legend(loc="lower right");

    #
    #
    # Non-linear response to Vcmax
    #
    #
    _canopy = CanopyLayer{FT}(n_leaf=1, g_max25=0.3);
    _envir  = AirLayer{FT}();
    _hs     = LeafHydraulics{FT}();
    _psm    = C3CLM(FT);
    _sm     = OSMWang{FT}();
    _canopy.APAR .= 500;
    _canopy.p_ups = -1.0;
    update_leaf_TP!(_psm, _canopy, _hs, _envir);
    temperature_effects!(_hs, _canopy.T);

    _vcmaxs = collect(FT, 10:2:100);
    _ans    = similar(_vcmaxs);
    _gsws   = similar(_vcmaxs);
    @showprogress for _i in eachindex(_vcmaxs)
        _canopy.APAR      .= _vcmaxs[_i] * 5;
        _canopy.ps.Vcmax25 = _vcmaxs[_i];
        _canopy.ps.Jmax25  = _vcmaxs[_i] * 1.67;
        _canopy.ps.Rd25    = _vcmaxs[_i] * 0.015;
        _canopy.ps.Vcmax   = _vcmaxs[_i];
        _canopy.ps.Jmax    = _vcmaxs[_i] * 1.67;
        _canopy.ps.Rd      = _vcmaxs[_i] * 0.015;
        _canopy.ps.APAR    = _apars[_i];
        leaf_ETR!(_psm, _canopy.ps);
        for _n in 1:600
            gas_exchange!(_psm, _canopy, _envir, GswDrive());
            for _iLF in eachindex(_canopy.An)
                _a1 = _canopy.Ac[_iLF];
                _a2 = _canopy.Aj[_iLF];
                _a3 = _canopy.Ap[_iLF];
                _ag = lower_quadratic(0.99, -(_a1 + _a2), _a1 * _a2);
                _ag = lower_quadratic(0.99, -(_ag + _a3), _ag * _a3);
                _canopy.Ag[_iLF] = _ag;
                _canopy.An[_iLF] = _ag - _canopy.ps.Rd;
            end;
            update_gsw!(_canopy, _sm, _psm, _envir, FT(60); τ=FT(1e-6), smoothing=true);
            gsw_control!(_psm, _canopy, _envir);
        end;
        #gas_exchange!(_psm, _canopy, _hs, _envir, _sm);
        for _iLF in eachindex((_canopy).An)
            _a1 = (_canopy).Ac[_iLF];
            _a2 = (_canopy).Aj[_iLF];
            _a3 = (_canopy).Ap[_iLF];
            _ag = lower_quadratic(0.99, -(_a1 + _a2), _a1 * _a2);
            _ag = lower_quadratic(0.99, -(_ag + _a3), _ag * _a3);
            (_canopy).Ag[_iLF] = _ag;
            (_canopy).An[_iLF] = _ag - (_canopy).ps.Rd;
        end;
        _ans[_i]  = _canopy.An[1];
        _gsws[_i] = _canopy.g_sw[1];
        _ans[_i]  = _canopy.An[1];
        _gsws[_i] = _canopy.g_sw[1];
    end

    # save figure
    _tx2 = _ax2.twinx();
    _ax2.plot(_vcmaxs, _ans, "g-", label=LS_Anet);
    _ax2.plot(0, 0, "c-", label="\$g_\\text{sw}\$");
    _tx2.plot(_vcmaxs, _gsws, "c-");

    _ax2.plot([14,14], [-1,22], "k:");
    _ax2.plot([78,78], [-1,22], "k:");

    _ax2.plot([14,46,78], [_ans[3],(_ans[3]+_ans[35])/2,_ans[35]], "g:");
    _ax2.plot([14,46,78], [_ans[3],(_ans[3]+_ans[35])/2,_ans[35]], "go");
    _ax2.plot(46, _ans[19], "go", mfc="none");

    _tx2.plot([14,46,78], [_gsws[3],(_gsws[3]+_gsws[35])/2,_gsws[35]], "c:");
    _tx2.plot([14,46,78], [_gsws[3],(_gsws[3]+_gsws[35])/2,_gsws[35]], "co");
    _tx2.plot(46, _gsws[19], "co", mfc="none");

    _ax2.annotate("", xy=(46,_ans[19]*0.98) , xytext=(46,(_ans[3]+_ans[35])/2*1.02)  , arrowprops=ARROW_G, fontsize=16, color="g");
    _tx2.annotate("", xy=(46,_gsws[19]*0.98), xytext=(46,(_gsws[3]+_gsws[35])/2*1.02), arrowprops=ARROW_C, fontsize=16, color="c");

    set_titles!(_axs; loc="left");
    set_xlabels!(_axs, [LS_APAR_unit, LS_Vcmax_unit]);
    set_ylabels!([_ax1,_tx1], [LS_Anet_unit, LS_gsw_unit]);
    set_ylabels!([_ax2,_tx2], [LS_Anet_unit, LS_gsw_unit]);
    set_ylims!([_ax1,_tx1,_ax2,_tx2], [[-1,22], [0,0.32], [-1,22], [0,0.32]]);
    save_canvas!(_fig, "figures/2021_canopy_complexity/2-nl.pdf", saving);
    if saving close(_fig) end;

    return nothing
end




"""
    plot_CC2021_APAR!(
                proj::CanopyComplexity2021{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot APAR comparison among canopy complexity levels, given
- `proj` [`CanopyComplexity2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_CC2021_APAR!(
            proj::CanopyComplexity2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    # use latex and serif
    if use_latex use_serif_tex(); end

    # create nodes and run simulations
    _nodes, _envir = create_spac(proj, OSMWang{FT}(), true);
    simulation!(proj, _nodes, _envir; τ=FT(1e-6));

    # create canvas and plot results
    _fig,_axs = create_canvas("3-apar.pdf");
    _ax1, = _axs;

    _ijkx_max = [maximum(_nodes[1].plant_ps[_i].APAR[1:end-1]) for _i in eachindex(_nodes[1].plant_ps)];
    _ijkx_min = [minimum(_nodes[1].plant_ps[_i].APAR[1:end-1]) for _i in eachindex(_nodes[1].plant_ps)];
    _2kx_sl   = [_nodes[2].plant_ps[_i].APAR[1] for _i in eachindex(_nodes[2].plant_ps)];
    _2kx_sh   = [_nodes[2].plant_ps[_i].APAR[2] for _i in eachindex(_nodes[2].plant_ps)];
    _kx       = [_nodes[3].plant_ps[_i].APAR[1] for _i in eachindex(_nodes[3].plant_ps)];
    _2x_sl    = [_nodes[4].plant_ps[1].APAR[1]];
    _2x_sh    = [_nodes[4].plant_ps[2].APAR[1]];
    _1x       = [_nodes[5].plant_ps[1].APAR[1]];
    _ax1.plot(_1x, [10.5], "o", color=COLORS[5], label="1X");
    _ax1.plot(_2x_sh, [10.5], "o", color=COLORS[4], mfc="none", label="2X sh");
    _ax1.plot(_2x_sl, [10.5], "o", color=COLORS[4], label="2X sl");
    _ax1.plot(_kx, eachindex(_kx), "-", color=COLORS[3], label="KX");
    _ax1.plot(_2kx_sl, eachindex(_2kx_sl), "-", color=COLORS[2], label="2KX sl");
    _ax1.plot(_2kx_sh, eachindex(_2kx_sl), ":", color=COLORS[2], label="2KX/IJKX sh");
    _ax1.fill_betweenx(eachindex(_2kx_sl), _ijkx_min, _ijkx_max, color=COLORS[1], label="IJKX sl", alpha=0.3);
    _ax1.legend(loc="lower right");

    set_xylims!(_axs, [0,720], [0,21]);
    set_xyticks!(_axs, [0,100,200,300,400,500,600,700], [0,4,8,12,16,20]);
    set_xylabels!(_axs, LS_APAR_unit, "Layer from bottom");
    save_canvas!(_fig, "figures/2021_canopy_complexity/3-apar.pdf", saving);
    if saving close(_fig) end;

    return nothing
end




"""
    plot_CC2021_SA!(
                proj::CanopyComplexity2021{FT};
                saving::Bool = false,
                use_latex::Bool = true,
                gradients::Bool = true
    ) where {FT<:AbstractFloat}

Plot canopy complexity sensitivity analysis, given
- `proj` [`CanopyComplexity2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
- `gradients` Optional. If true, apply a vertical canopy gradient
"""
function plot_CC2021_SA!(
            proj::CanopyComplexity2021{FT};
            saving::Bool = false,
            use_latex::Bool = true,
            gradients::Bool = false
) where {FT<:AbstractFloat}
    # use latex and serif
    if use_latex use_serif_tex(); end

    # create nodes and envir to work on
    _nodes, _envir = create_spac(proj, OSMWang{FT}(), gradients);
    _legends = ["IJKX", "2KX", "KX", "2X", "1X"];
    _alphas  = [0.7, 0.7, 0.7, 0.7, 0.7];
    _lines   = ["--", "-", "-", "-", "-"];

    # use the same mean Vcmax
    if !gradients
        for _node in _nodes
            update_VJRWW!(_node, FT(35.30443317200483));
        end
    end

    # a figure to plot all simulations
    @info tinfo("Plot sensitivity analysis...");
    _id = gradients ? "CC-4" : "CC-3";
    _fn = gradients ? "5-sa.pdf" : "4-sa.pdf";
    _fig,_axs = create_canvas(_id; nrow=3, ncol=5, figsize=(12,9));

    # a figure to plot all simulations in three panels
    _error_2kx_CO₂ = FT[];
    _error_2kx_H₂O = FT[];
    _error_2kx_SIF = FT[];
    _error_2x_CO₂  = FT[];
    _error_2x_H₂O  = FT[];
    _error_2x_SIF  = FT[];

    #
    #
    # models responses to radiation
    #
    #
    _trad    = numerical∫(_nodes[1].in_rad.E_direct , _nodes[1].wl_set.dWL) + numerical∫(_nodes[1].in_rad.E_diffuse, _nodes[1].wl_set.dWL);
    _rrads   = collect(FT, 0.1:0.1:1);
    _mat_CO₂ = zeros(FT, (length(_rrads),5));
    _mat_H₂O = zeros(FT, (length(_rrads),5));
    _mat_SIF = zeros(FT, (length(_rrads),5));
    @showprogress for _i in eachindex(_rrads)
        _env_tmp = deepcopy(_envir);
        _nds_tmp = deepcopy(_nodes);
        for _node in _nds_tmp
            _node.in_rad.E_direct  .*= _rrads[_i];
            _node.in_rad.E_diffuse .*= _rrads[_i];
        end
        simulation!(proj, _nds_tmp, _env_tmp; τ=FT(1e-6));
        _mat_CO₂[_i,:] .= [_nds_tmp[_j].f_npp for _j in eachindex(_nds_tmp)];
        _mat_H₂O[_i,:] .= [_nds_tmp[_j].f_H₂O for _j in eachindex(_nds_tmp)];
        _mat_SIF[_i,:] .= [SIF_740(_nds_tmp[_j].can_rad, _nds_tmp[_j].wl_set) for _j in eachindex(_nds_tmp)];
    end
    @inline plot_subplots(ax1,ax2,ax3,xs) = (
        for _mode in eachindex(_legends)
            ax1.plot(xs, _mat_CO₂[:,_mode]        , label=_legends[_mode], alpha=_alphas[_mode], linestyle=_lines[_mode]);
            ax2.plot(xs, _mat_H₂O[:,_mode] .* 1000, label=_legends[_mode], alpha=_alphas[_mode], linestyle=_lines[_mode]);
            ax3.plot(xs, _mat_SIF[:,_mode]        , label=_legends[_mode], alpha=_alphas[_mode], linestyle=_lines[_mode]);
        end;
        push!(_error_2kx_CO₂, (_mat_CO₂[:,2] ./ _mat_CO₂[:,1])...);
        push!(_error_2kx_H₂O, (_mat_H₂O[:,2] ./ _mat_H₂O[:,1])...);
        push!(_error_2kx_SIF, (_mat_SIF[:,2] ./ _mat_SIF[:,1])...);
        push!(_error_2x_CO₂ , (_mat_CO₂[:,4] ./ _mat_CO₂[:,1])...);
        push!(_error_2x_H₂O , (_mat_H₂O[:,4] ./ _mat_H₂O[:,1])...);
        push!(_error_2x_SIF , (_mat_SIF[:,4] ./ _mat_SIF[:,1])...);
    );
    plot_subplots(_axs[1], _axs[6], _axs[11], _rrads .* _trad ./ 1000);

    #
    #
    # models responses to temperature
    #
    #
    _tems    = collect(FT, 283.1:3:313.1);
    _mat_CO₂ = zeros(FT, (length(_tems),5));
    _mat_H₂O = zeros(FT, (length(_tems),5));
    _mat_SIF = zeros(FT, (length(_tems),5));
    @showprogress for _i in eachindex(_tems)
        _env_tmp       = deepcopy(_envir);
        _env_tmp.t_air = _tems[_i];
        _env_tmp.p_sat = saturation_vapor_pressure(_tems[_i]);
        _env_tmp.p_H₂O = _env_tmp.p_sat * _env_tmp.RH;
        _env_tmp.vpd   = _env_tmp.p_sat - _env_tmp.p_H₂O;
        _nds_tmp       = deepcopy(_nodes);
        for _node in _nds_tmp
            for _j in eachindex(_node.plant_ps)
                _iHS   = _node.plant_hs.leaves[_j];
                _iPS   = _node.plant_ps[_j];
                _iPS.T = _tems[_i];
                update_leaf_TP!(_node.photo_set, _iPS, _iHS, _env_tmp);
                temperature_effects!(_iHS, _tems[_i]);
            end
        end
        simulation!(proj, _nds_tmp, _env_tmp; τ=FT(1e-6));
        _mat_CO₂[_i,:] .= [_nds_tmp[_j].f_npp for _j in eachindex(_nds_tmp)];
        _mat_H₂O[_i,:] .= [_nds_tmp[_j].f_H₂O for _j in eachindex(_nds_tmp)];
        _mat_SIF[_i,:] .= [SIF_740(_nds_tmp[_j].can_rad, _nds_tmp[_j].wl_set) for _j in eachindex(_nds_tmp)];
    end
    plot_subplots(_axs[2], _axs[7], _axs[12], _tems);

    #
    #
    # models responses to VPD
    #
    #
    _vpds    = collect(FT, 500:350:4000);
    _mat_CO₂ = zeros(FT, (length(_vpds),5));
    _mat_H₂O = zeros(FT, (length(_vpds),5));
    _mat_SIF = zeros(FT, (length(_vpds),5));
    @showprogress for _i in eachindex(_vpds)
        _env_tmp       = deepcopy(_envir);
        _env_tmp.vpd   = _vpds[_i];
        _env_tmp.p_H₂O = _env_tmp.p_sat - _vpds[_i];
        _env_tmp.RH    = _env_tmp.p_H₂O / _env_tmp.p_sat;
        _nds_tmp       = deepcopy(_nodes);
        simulation!(proj, _nds_tmp, _env_tmp; τ=FT(1e-6));
        _mat_CO₂[_i,:] .= [_nds_tmp[_j].f_npp for _j in eachindex(_nds_tmp)];
        _mat_H₂O[_i,:] .= [_nds_tmp[_j].f_H₂O for _j in eachindex(_nds_tmp)];
        _mat_SIF[_i,:] .= [SIF_740(_nds_tmp[_j].can_rad, _nds_tmp[_j].wl_set) for _j in eachindex(_nds_tmp)];
    end
    plot_subplots(_axs[3], _axs[8], _axs[13], _vpds);

    #
    #
    # models responses to Psoil
    #
    #
    _psoils  = collect(FT, -0:-0.5:-5.0);
    _mat_CO₂ = zeros(FT, (length(_psoils),5));
    _mat_H₂O = zeros(FT, (length(_psoils),5));
    _mat_SIF = zeros(FT, (length(_psoils),5));
    @showprogress for _i in eachindex(_psoils)
        _env_tmp = deepcopy(_envir);
        _nds_tmp = deepcopy(_nodes);
        for _node in _nds_tmp
            for _root in _node.plant_hs.roots
                _root.p_ups = _psoils[_i];
            end
        end
        simulation!(proj, _nds_tmp, _env_tmp; τ=FT(1e-6));
        _mat_CO₂[_i,:] .= [_nds_tmp[_j].f_npp for _j in eachindex(_nds_tmp)];
        _mat_H₂O[_i,:] .= [_nds_tmp[_j].f_H₂O for _j in eachindex(_nds_tmp)];
        _mat_SIF[_i,:] .= [SIF_740(_nds_tmp[_j].can_rad, _nds_tmp[_j].wl_set) for _j in eachindex(_nds_tmp)];
    end
    plot_subplots(_axs[4], _axs[9], _axs[14], _psoils);

    #
    #
    # models responses to CO₂
    #
    #
    _pas     = collect(FT, 20:5:80);
    _mat_CO₂ = zeros(FT, (length(_pas),5));
    _mat_H₂O = zeros(FT, (length(_pas),5));
    _mat_SIF = zeros(FT, (length(_pas),5));
    @showprogress for _i in eachindex(_pas)
        _env_tmp     = deepcopy(_envir);
        _env_tmp.p_a = _pas[_i];
        _nds_tmp     = deepcopy(_nodes);
        simulation!(proj, _nds_tmp, _env_tmp; τ=FT(1e-6));
        _mat_CO₂[_i,:] .= [_nds_tmp[_j].f_npp for _j in eachindex(_nds_tmp)];
        _mat_H₂O[_i,:] .= [_nds_tmp[_j].f_H₂O for _j in eachindex(_nds_tmp)];
        _mat_SIF[_i,:] .= [SIF_740(_nds_tmp[_j].can_rad, _nds_tmp[_j].wl_set) for _j in eachindex(_nds_tmp)];
    end
    plot_subplots(_axs[5], _axs[10], _axs[15], _pas);

    _axs[1].legend(loc="lower right");
    set_xlabels!(_axs[11:15], ["Radiation (W m\$^{-2}\$)", "\$T\$ (K)", "VPD (Pa)", "\$\\upPsi_\\text{soil}\$ (MPa)", "\$P_\\text{CO2}\$ (Pa)"]);
    set_ylabels!(_axs[[1,6,11]], [LS_NEE_unit, LS_ET_unit_m, LS_SIF740_unit]);
    set_titles!(_axs[1:5]; loc="left");
    save_canvas!(_fig, "figures/2021_canopy_complexity/$(_fn)", saving);
    if saving close(_fig) end;

    # display error information
    @show nanmean(_error_2kx_CO₂), nanmean(_error_2x_CO₂);
    @show nanmean(_error_2kx_H₂O), nanmean(_error_2x_H₂O);
    @show nanmean(_error_2kx_SIF), nanmean(_error_2x_SIF);
    @show nanstd(_error_2kx_CO₂), nanstd(_error_2x_CO₂);
    @show nanstd(_error_2kx_H₂O), nanstd(_error_2x_H₂O);
    @show nanstd(_error_2kx_SIF), nanstd(_error_2x_SIF);

    return nothing
end




"""
    plot_CC2021_SA_err!(
                proj::CanopyComplexity2021{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot CanopyComplexity2021 error comparison figure, given
- `proj` [`CanopyComplexity2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_CC2021_SA_err!(
            proj::CanopyComplexity2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    # use latex and serif
    if use_latex use_serif_tex(); end

    # hard code the errors
    _locs = [1,2,3,4,6,7,8,9,11,12,13,14];
    _errs = [2.36,11.1,5.44,23.4,1.19,3.67,3.76,8.16,2.78,7.87,4.16,13.3];
    _stds = [1.15,3.77,3.39,9.17,0.89,3.07,3.62,7.96,1.57,5.56,5.33,8.44];
    _cols = repeat(["r", "r", "c", "c"], 3);

    # plot the figure
    @info tinfo("Plot error comparison...");
    _fig,_axs = create_canvas("CC-7"; figsize=(5,3.5));
    _ax1 = _axs[1];
    _i1  = [1,5, 9];
    _i2  = [2,6,10];
    _i3  = [3,7,11];
    _i4  = [4,8,12];
    _ax1.bar(_locs[_i1], _errs[_i1], yerr=_stds[_i1], color=_cols[_i1], alpha=0.5, label="2KX + uniform \$V_\\text{cmax}\$"         );
    _ax1.bar(_locs[_i2], _errs[_i2], yerr=_stds[_i2], color=_cols[_i2], alpha=0.9, label="2KX + vertical \$V_\\text{cmax}\$ profile");
    _ax1.bar(_locs[_i3], _errs[_i3], yerr=_stds[_i3], color=_cols[_i3], alpha=0.5, label="2X + uniform \$V_\\text{cmax}\$"          );
    _ax1.bar(_locs[_i4], _errs[_i4], yerr=_stds[_i4], color=_cols[_i4], alpha=0.9, label="2X + vertical \$V_\\text{cmax}\$ profile" );
    _ax1.legend(loc="upper right");
    set_xticks!(_axs, [2.5,7.5,12.5]);
    set_xticklabels!(_axs, ["NEE", "ET", "SIF\$_\\text{740}\$"]);
    set_ylabels!(_axs, "Relative difference (\\%)");
    save_canvas!(_fig, "figures/2021_canopy_complexity/6-err.pdf", saving);
    if saving close(_fig) end;

    return nothing
end




"""
    plot_CC2021_SIF_RAD!(
                proj::CanopyComplexity2021{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot CanopyComplexity2021 SIF vs. radiation discussion figure, given
- `proj` [`CanopyComplexity2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_CC2021_SIF_RAD!(
            proj::CanopyComplexity2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    # use latex and serif
    if use_latex use_serif_tex(); end

    # PSII yield and SIF yield at leaf level
    _nodes,  = create_spac(proj, OSMWang{FT}(), true);
    _leaf    = Leaf{FT}();
    _psm     = C3CLM(FT);
    _envir   = AirLayer{FT}();
    _pis     = collect(FT, 5:1:60);
    _apars   = collect(FT, 10:10:500);
    _vcmaxs  = [_ps.ps.Vcmax25 for _ps in _nodes[3].plant_ps];
    _mat_psy = zeros(FT, length(_pis), length(_apars));
    _mat_fqe = zeros(FT, length(_pis), length(_apars));
    _mat_vcm = zeros(FT, length(_vcmaxs), length(_apars));
    for _i in eachindex(_pis), _j in eachindex(_apars)
        _leaf.APAR = _apars[_j];
        leaf_photosynthesis!(_psm, _leaf, _envir, PCO₂Mode(), _pis[_i]);
        leaf_fluorescence!(_psm.Flu, _leaf);
        _mat_psy[_i,_j] = _leaf.φ;
        _mat_fqe[_i,_j] = _leaf.φs * 100;
    end;
    for _i in eachindex(_vcmaxs), _j in eachindex(_apars)
        _leaf.p_i     = 20;
        _leaf.APAR    = _apars[_j];
        _leaf.Vcmax25 = _vcmaxs[_i];
        _leaf.Jmax25  = _vcmaxs[_i] * 1.8;
        _leaf.Rd25    = _vcmaxs[_i] * 0.015;
        _leaf.T_old   = 0;
        leaf_photosynthesis!(_psm, _leaf, _envir, PCO₂Mode(), _pis[_i]);
        leaf_fluorescence!(_psm.Flu, _leaf);
        _mat_vcm[_i,_j] = _leaf.φs * 100;
    end;

    # plot the results
    _fig,_axs = create_canvas("CC-PSY-FQE"; nrow=3, ncol=2, figsize=(7.8,9.5));
    _ax1,_ax2,_ax3,_ax4,_ax5,_ax6 = _axs;
    _cm1 = _ax1.pcolor(_pis, _apars, _mat_psy'; shading="auto");
    _cm2 = _ax2.pcolor(_pis, _apars, _mat_fqe'; shading="auto");
    for _i in [6,16,26,36]
        _ax3.plot(_apars, _apars .* _mat_fqe[_i,:] ./ 100);
        _ax3.text(_apars[end], _apars[end] * _mat_fqe[_i,end] / 100, "$(_pis[_i])", ha="left", va="center");
    end;
    for _i in [5,15,25,30,35,45]
        _ax4.plot(_pis, _apars[_i] .* _mat_fqe[:,_i] ./ 100, label="APAR = $(_apars[_i])");
        _ax4.text(_pis[end], _apars[_i] * _mat_fqe[end,_i] / 100, "$(_apars[_i])", ha="right", va="bottom");
    end;

    # plot the complicated scenario without Vcmax profile
    _nodes, = create_spac(proj, OSMWang{FT}(), false);
    for _pi in [20,30,40,50]
        @info tinfo("Simulating the case at a CO₂ of $(_pi) Pa...");
        _env_tmp     = deepcopy(_envir);
        _env_tmp.p_a = _pi;
        _nds_tmp     = deepcopy(_nodes);
        simulation!(proj, _nds_tmp, _env_tmp; τ=FT(1e-6));
        _sim_fqe = [_ps.φs[1] * 100 for _ps in _nds_tmp[3].plant_ps];
        _ax5.plot(_sim_fqe, eachindex(_vcmaxs), label="$(_pi) Pa");
    end;
    _ax5.legend(loc="lower right");

    # plot the complicated scenario with Vcmax profile
    _nodes, = create_spac(proj, OSMWang{FT}(), true);
    for _pi in [20,30,40,50]
        @info tinfo("Simulating the case at a CO₂ of $(_pi) Pa...");
        _env_tmp     = deepcopy(_envir);
        _env_tmp.p_a = _pi;
        _nds_tmp     = deepcopy(_nodes);
        simulation!(proj, _nds_tmp, _env_tmp; τ=FT(1e-6));
        _sim_fqe = [_ps.φs[1] * 100 for _ps in _nds_tmp[3].plant_ps];
        _ax6.plot(_sim_fqe, eachindex(_vcmaxs), label="$(_pi) Pa");
    end;
    _ax6.legend(loc="upper left");

    _fig.colorbar(_cm1; ax=_ax1, label="PSII yield (-)");
    _fig.colorbar(_cm2; ax=_ax2, label="\$\\phi_\\text{F}\$ (\\%)");
    set_titles!(_axs; loc="left");
    set_xlabels!(_axs, ["Internal CO\$_2\$ (Pa)"; "Internal CO\$_2\$ (Pa)"; LS_APAR_unit; "Internal CO\$_2\$ (Pa)"; "\$\\phi_\\text{F}\$ (\\%)"; "\$\\phi_\\text{F}\$ (\\%)"]);
    set_ylabels!(_axs, [LS_APAR_unit; LS_APAR_unit; repeat(["\$\\phi_\\text{F}\$ \$\\cdot\$ $(LS_APAR_unit)"],2); repeat(["Layer from bottom"],2)]);
    save_canvas!(_fig, "figures/2021_canopy_complexity/11-psy-fqe.pdf", saving);
    if saving close(_fig) end;

    return nothing
end
