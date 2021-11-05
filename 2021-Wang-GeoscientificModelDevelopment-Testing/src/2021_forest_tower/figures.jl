###############################################################################
#
# Fig. 2: Spectrum and APAR differences for ForestTower2021 project
#
###############################################################################
"""
    plot_FT2021_spectrum!(
                proj::ForestTower2021{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot spectrum and APAR differences for project ForestTower2021, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_FT2021_spectrum!(
            proj::ForestTower2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    @info tinfo("Plot difference in spectrum and APAR...");

    # use latex and serif
    if use_latex use_serif_tex(); end

    # run the canopy RT at default setup
    collections = initialize_rt_module(FT; nLayer=20, LAI=3);
    angles, can, can_opt, can_rad, in_rad, leaves, rt_con, rt_dim, soil,
            wls = collections;
    can_rad_car = deepcopy(can_rad);
    leaf_car    = deepcopy(leaves[1]);

    # run the model again without Car absorption for APAR
    for leaf in leaves
        fluspect!(leaf, wls; APAR_car=false);
    end
    canopy_geometry!(can, angles, can_opt, rt_con);
    canopy_matrices!(leaves, can_opt);
    short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
    canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
    can_rad_noc = deepcopy(can_rad);
    leaf_noc    = deepcopy(leaves[1]);

    # calculate the APAR difference matrix for the first layer
    apar_sunlit_diff = (can_rad_car.absPAR_sunCab[:,:,1] .-
                        can_rad_noc.absPAR_sunCab[:,:,1]) .* 1e6;
    apar_shaded_diff = (can_rad_car.absPAR_shadeCab[1] .-
                        can_rad_noc.absPAR_shadeCab[1]) .* 1e6;
    apar_sunlit_diff = [apar_sunlit_diff apar_sunlit_diff[:,1,1]];

    # create canvas
    fig = figure("FT-spectrum"; figsize=(6.5,3.5), dpi=100);
    fig.clear();
    ax1 = fig.add_subplot(1, 2, 1);
    ax2 = fig.add_subplot(1, 2, 2, polar=true);

    ax1.plot(wls.WL, leaf_car.kChlrel, "k-", label="Chl+Car" );
    ax1.plot(wls.WL, leaf_noc.kChlrel, "k:", label="Chl only");
    ax1.legend(loc="lower left");

    cm = ax2.contourf(deg2rad.(collect(5:10:370)), collect(5:10:85),
                      apar_sunlit_diff);
    cb = fig.colorbar(cm, ax=ax2, label="\$\\upDelta\$A" * LS_PAR_unit);
    cb.ax.plot([0,400], [apar_shaded_diff for i in 1:2], "w-");
    ax2.set_yticks([]);
    ax2.grid(linestyle=":");

    set_xlims!([ax1], [390,810]);
    set_ylims!([ax1,ax2], [[0,1], [0,85]]);
    set_titles!([ax1,ax2]; loc="left");
    set_xylabels!([ax1], "Wavelength (nm)", "Relative absorption (-)");

    save_canvas!(fig, "figures/2021_forest_tower/2_spectrum.pdf", saving);
    if saving close(fig) end;

    return nothing;
end








###############################################################################
#
# Fig. 3: Impact of SZA and LAI on the matricies above
#
###############################################################################
"""
    plot_FT2021_clumping!(
                proj::ForestTower2021{FT};
                saving::Bool=false,
                use_latex::Bool=true
    ) where {FT<:AbstractFloat}

Plot clumping index impacts on canopy radiative transfter, given
- `proj` [`ForestTower2021`](@ref) type project control
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font

Parameterization of the simulation
- Canopy layers = 30
- Solar zenith angle = 30
- Leaf area index = 3.0
- Clumping index in [0.5, 1.0]
- VZA = 0, RAA = 0 (defaults)

Key components of the results
- Sunlit leaf fraction at each canopy layer
- Mean sunlit leaf APAR for each canopy layer
- Mean shaded leaf APAR for each canopy layer
"""
function plot_FT2021_clumping!(
            proj::ForestTower2021{FT};
            saving::Bool=false,
            use_latex::Bool=true
) where {FT<:AbstractFloat}
    # use serif and latex wrapper
    if use_latex use_serif_tex(); end

    @info tinfo("Plotting clumping index impacts...");

    # create canvas to plot the results
    _fig,_axs = create_canvas("FT-clumping"; ncol=2, nrow=2);
    _ax1,_ax2,_ax3,_ax4 = _axs;

    # figure attributes
    _cis     = [1.0, 0.5];
    _leg_ci  = ["CI = $(_ci)" for _ci in _cis];
    _leg_csl = ["CI = $(_ci) sunlit" for _ci in _cis];
    _leg_csh = ["CI = $(_ci) shaded" for _ci in _cis];
    _lines   = ["k:", "k-"];
    _lines2  = ["k--", "k-"];
    _xlabels = ["Sunlit leaf fraction (-)", LS_APAR_unit, "Wavelength (nm)",
                "Wavelength (nm)"];
    _ylabels = ["\$i\$th canopy layer", "\$i\$th canopy layer", "Albedo (-)",
                "SIF " * LS_unit_SIF]

    # plot the sensitivity to CI
    _angles, _can, _can_opt, _can_rad, _in_rad, _leaves, _rt_con, _rt_dim,
             _soil, _wls = initialize_rt_module(FT; nLayer=30, LAI=FT(3));
    for _i in eachindex(_cis)
        _can.Ω = _cis[_i];
        canopy_geometry!(_can, _angles, _can_opt, _rt_con);
        canopy_matrices!(_leaves, _can_opt);
        short_wave!(_can, _can_opt, _can_rad, _in_rad, _soil, _rt_con);
        canopy_fluxes!(_can, _can_opt, _can_rad, _in_rad, _soil, _leaves, _wls,
                       _rt_con);
        SIF_fluxes!(_leaves, _can_opt, _can_rad, _can, _soil, _wls, _rt_con,
                    _rt_dim);

        # calculate the absorbed solar radiation
        _fs  = (_can_opt.Ps[1:end-1] .+ _can_opt.Ps[2:end]) ./ 2;
        _esh = _can_rad.intNetSW_shade  .* (1 .- _fs) .* 0.1;
        _esl = _can_rad.intNetSW_sunlit .* _fs .* 0.1;
        _ess = _esl + _esh;

        # plot the model simulations
        _ax1.plot(_can_opt.Ps, collect(31:-1:1) .- 0.5, _lines[_i],
                  label=_leg_ci[_i]);
        _ax2.plot([nanmean(_can_rad.absPAR_sunCab[:,:,_l])*1e6 for _l in 1:30],
                  collect(30:-1:1), _lines[_i], label=_leg_csl[_i]);
        _ax2.plot(_can_rad.absPAR_shadeCab .* 1e6, collect(30:-1:1),
                  _lines2[_i], label=_leg_csh[_i], alpha=0.5);
        _ax3.plot(_wls.WL, _can_rad.alb_obs, _lines[_i], label=_leg_ci[_i]);
        _ax4.plot(_wls.WLF, _can_rad.SIF_obs, _lines[_i], label=_leg_ci[_i]);
    end

    # set up axis attributes
    _ax1.legend(loc="lower right");
    _ax2.legend(loc="lower right");
    _ax3.legend(loc="upper right");
    _ax4.legend(loc="upper right");
    set_xylabels!([_ax1,_ax2,_ax3,_ax4], _xlabels, _ylabels);
    set_titles!(_axs; loc="left");

    # save the figure
    save_canvas!(_fig, "figures/2021_forest_tower/3_clumping.pdf", saving);
    if saving close(_fig) end;

    return nothing
end








###############################################################################
#
# Fig. 4: Supply curve and Ecrit
#
###############################################################################
"""
    plot_FT2021_ecrit!(
                proj::ForestTower2021{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot Ecrit examples, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_FT2021_ecrit!(
            proj::ForestTower2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    @info tinfo("Plot Ecrit calculations...");

    # use latex and serif
    if use_latex use_serif_tex(); end

    # create a figure
    _fig,_axs = create_canvas("FT-Ecrit", ncol=2);
    _ax1,_ax2 = _axs;

    # create a supply curve with p_ups = 0
    _leaf = LeafHydraulics{FT}(N=10000);
    _e  = FT(0);
    _es = FT[];
    _ps = FT[];
    while true
        _p = end_pressure(_leaf, _e);
        push!(_es, _e);
        push!(_ps, _p);
        if _p < -4
            break
        end
        _e += FT(1e-4);
    end
    _ax1.plot(_ps, _es, "k-");

    # create a supply curve with P_ups = -1
    _leaf.p_ups = -1;
    _e  = FT(0);
    _es = FT[];
    _ps = FT[];
    while true
        _p = end_pressure(_leaf, _e);
        push!(_es, _e);
        push!(_ps, _p);
        if _p < -4
            break
        end
        _e += FT(1e-4);
    end
    _ax1.plot(_ps, _es, "k--");
    _ax1.plot([_leaf.p_crt, _leaf.p_crt], [-0.001, 0.075], "k-", alpha=0.5);

    # plot e_crit as a function of p_ups
    _p_upss = collect(FT, 0:-0.1:-3);
    _e_crts = FT[];
    for _p_ups in _p_upss
        _leaf.p_ups = _p_ups;
        _e_crt = critical_flow(_leaf);
        push!(_e_crts, _e_crt);
    end
    _ax2.plot(_p_upss, _e_crts, "k-");

    # set the axis attributes
    set_xylims!(_axs, [[-3.5,0], [-3,0]], [-0.001,0.075]);
    set_xylabels!(_axs, "\$\\upPsi\$ (MPa)", [LS_E_unit, LS_Ecrit_unit]);
    set_titles!(_axs; loc="left");

    save_canvas!(_fig, "figures/2021_forest_tower/4_ecrit.pdf", saving);
    if saving close(_fig) end;

    return nothing
end








###############################################################################
#
# Figs. 5: Examples of model responses to the environment
#
###############################################################################
"""
    plot_FT2021_responses!(
                proj::ForestTower2021{FT};
                saving::Bool = false,
                use_latex::Bool = true,
                site::String = "NiwotRidge"
    ) where {FT<:AbstractFloat}
Plot stomatal conductance responses to different environmental stimuli, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_FT2021_responses!(
            proj::ForestTower2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    @info tinfo("Plot main text model responses figure...");

    # use latex and serif
    if use_latex use_serif_tex(); end

    # create canvas to plot on
    _fig,_axs = create_canvas("FT-responses", nrow=2, ncol=2);
    _gx1,_gx2,_cx1,_cx2 = _axs;

    # 1.1 test stomatal responses to VPD from 600 Pa to 2500 Pa
    _node_bbb = create_spac(proj, "NiwotRidge", ESMBallBerry{FT}());
    _node_bbm = create_spac(proj, "NiwotRidge", ESMBallBerry{FT}());
    _node_meb = create_spac(proj, "NiwotRidge", ESMMedlyn{FT}());
    _node_med = create_spac(proj, "NiwotRidge", ESMMedlyn{FT}());
    _node_osm = create_spac(proj, "NiwotRidge", OSMWang{FT}());

    # 1.2 generate data to plot
    _p_vpd = collect(FT, 500:100:2500);
    _c_bbb = similar(_p_vpd);
    _c_bbm = similar(_p_vpd);
    _c_meb = similar(_p_vpd);
    _c_med = similar(_p_vpd);
    _c_osm = similar(_p_vpd);
    _e_bbb = similar(_p_vpd);
    _e_bbm = similar(_p_vpd);
    _e_meb = similar(_p_vpd);
    _e_med = similar(_p_vpd);
    _e_osm = similar(_p_vpd);
    _ns    = [_node_bbm, _node_med, _node_osm, _node_bbb, _node_meb];
    _cs    = [_c_bbm, _c_med, _c_osm, _c_bbb, _c_meb];
    _es    = [_e_bbm, _e_med, _e_osm, _e_bbb, _e_meb];
    _bs    = [false, false, false, true, true];
    # initialize the nodes
    for _j in eachindex(_ns)
        _node = _ns[_j];
        for _envir in _node.envirs
            _envir.p_H₂O = _envir.p_sat - 500;
            _envir.vpd   = 500;
            _envir.RH    = _envir.p_H₂O / _envir.p_sat;
        end
        # make sure g_max is big enough so as to show the CO₂ response
        for _clayer in _node.plant_ps
            _clayer.g_max   = 0.2;
            _clayer.g_max25 = 0.2;
        end
        simulation!(proj, _node, _bs[_j]);
    end
    @showprogress for _i in eachindex(_p_vpd)
        _vpd = _p_vpd[_i];
        for _j in eachindex(_ns)
            _node = _ns[_j];
            _cmod = _cs[_j];
            _emod = _es[_j];
            for _envir in _node.envirs
                _envir.p_H₂O = _envir.p_sat - _vpd;
                _envir.vpd   = _vpd;
                _envir.RH    = _envir.p_H₂O / _envir.p_sat;
            end
            # make sure g_max is big enough so as to show the CO₂ response
            for _clayer in _node.plant_ps
                _clayer.g_max   = 0.2;
                _clayer.g_max25 = 0.2;
            end
            simulation!(proj, _node, _bs[_j]);
            _cmod[_i] = _node.f_npp / _node.ga;
            _emod[_i] = _node.f_H₂O / _node.ga * 1000;
        end
    end

    # 1.3 plot the simulation results
    _vpd = _p_vpd ./ _node_bbb.envirs[1].p_atm;
    _gx1.plot(_p_vpd, _e_bbm ./ _vpd, ":", color="tomato");
    _gx1.plot(_p_vpd, _e_med ./ _vpd, ":", color="royalblue");
    _gx1.plot(_p_vpd, _e_osm ./ _vpd, "-", color="c");
    _gx1.plot(_p_vpd, _e_bbb ./ _vpd, "-", color="tomato");
    _gx1.plot(_p_vpd, _e_meb ./ _vpd, "-", color="royalblue");
    _cx1.plot(_p_vpd, _c_bbm, ":", color="tomato");
    _cx1.plot(_p_vpd, _c_med, ":", color="royalblue");
    _cx1.plot(_p_vpd, _c_osm, "-", color="c");
    _cx1.plot(_p_vpd, _c_bbb, "-", color="tomato");
    _cx1.plot(_p_vpd, _c_meb, "-", color="royalblue");

    # 2.1 test stomatal responses to Psoil from 0 to -3 MPa
    _node_bbb = create_spac(proj, "NiwotRidge", ESMBallBerry{FT}());
    _node_bbm = create_spac(proj, "NiwotRidge", ESMBallBerry{FT}());
    _node_meb = create_spac(proj, "NiwotRidge", ESMMedlyn{FT}());
    _node_med = create_spac(proj, "NiwotRidge", ESMMedlyn{FT}());
    _node_osm = create_spac(proj, "NiwotRidge", OSMWang{FT}());

    # 2.2 generate data to plot
    _p_soil = collect(FT, 0.0:-0.2:-4.0);
    _c_bbb = similar(_p_soil);
    _c_bbm = similar(_p_soil);
    _c_meb = similar(_p_soil);
    _c_med = similar(_p_soil);
    _c_osm = similar(_p_soil);
    _e_bbb = similar(_p_soil);
    _e_bbm = similar(_p_soil);
    _e_meb = similar(_p_soil);
    _e_med = similar(_p_soil);
    _e_osm = similar(_p_soil);
    _ns    = [_node_bbm, _node_med, _node_osm, _node_bbb, _node_meb];
    _cs    = [_c_bbm, _c_med, _c_osm, _c_bbb, _c_meb];
    _es    = [_e_bbm, _e_med, _e_osm, _e_bbb, _e_meb];
    _bs    = [false, false, false, true, true];
    # initialize the nodes
    for _j in eachindex(_ns)
        _node = _ns[_j];
        for _root in _node.plant_hs.roots
            _root.p_ups = 0;
        end;
        # make sure g_max is big enough so as to show the CO₂ response
        for _clayer in _node.plant_ps
            _clayer.g_max   = 0.2;
            _clayer.g_max25 = 0.2;
        end
        simulation!(proj, _node, _bs[_j]);
    end
    @showprogress for _i in eachindex(_p_soil)
        _ps = _p_soil[_i];
        for _j in eachindex(_ns)
            _node = _ns[_j];
            _cmod = _cs[_j];
            _emod = _es[_j];
            # make sure g_max is big enough so as to show the CO₂ response
            for _clayer in _node.plant_ps
                _clayer.g_max   = 0.2;
                _clayer.g_max25 = 0.2;
            end
            for _root in _node.plant_hs.roots
                _root.p_ups = _ps;
            end;
            simulation!(proj, _node, _bs[_j]);
            _cmod[_i] = _node.f_npp / _node.ga;
            _emod[_i] = _node.f_H₂O / _node.ga * 1000;
        end
    end

    # 2.3 plot the simulation results
    _vpd = _node_bbb.envirs[1].vpd / _node_bbb.envirs[1].p_atm;
    _gx2.plot(_p_soil, _e_bbm ./ _vpd, ":", color="tomato");
    _gx2.plot(_p_soil, _e_med ./ _vpd, ":", color="royalblue");
    _gx2.plot(_p_soil, _e_osm ./ _vpd, "-", color="c");
    _gx2.plot(_p_soil, _e_bbb ./ _vpd, "-", color="tomato");
    _gx2.plot(_p_soil, _e_meb ./ _vpd, "-", color="royalblue");
    _cx2.plot(_p_soil, _c_bbm, ":", color="tomato", label="BBM no \$\\beta\$");
    _cx2.plot(_p_soil, _c_med, ":", color="royalblue", label="MED no \$\\beta\$");
    _cx2.plot(_p_soil, _c_osm, "-", color="c", label="OSM");
    _cx2.plot(_p_soil, _c_bbb, "-", color="tomato", label="BBM");
    _cx2.plot(_p_soil, _c_meb, "-", color="royalblue", label="MED");
    _cx2.legend(loc="lower right");

    # set axis attributes
    set_xlabels!([_cx1,_cx2], ["VPD (Pa)", "\$\\upPsi_\\text{soil}\$ (MPa)"]);
    set_ylabels!([_gx1,_cx1], ["\$G\$ (mmol m\$^{-2}\$ s\$^{-1}\$)", "CNPP" * latex_unit("A")]);
    set_titles!(_axs, loc="left");

    # save figure
    save_canvas!(_fig, "figures/2021_forest_tower/6_model_responses.pdf", saving);
    if saving close(_fig) end;

    return nothing
end








###############################################################################
#
# Figs. 7 and 8: Examples of model comparisons
#
###############################################################################
"""
    plot_FT2021_example!(
                proj::ForestTower2021{FT},
                nTH::Int;
                saving::Bool = false,
                use_latex::Bool = true,
                site::String = "NiwotRidge",
                g1::Bool = false
    ) where {FT<:AbstractFloat}

Plot model examples for project ForestTower2021 for a site, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `nTH` Figure number
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
- `site` Site of flux tower, must be `NiwotRidge` or `Ozark`
- `g1` If `true`, empirical g1 was fitted. Default is `false`
"""
function plot_FT2021_example!(
            proj::ForestTower2021{FT},
            nTH::Int;
            saving::Bool = false,
            use_latex::Bool = true,
            site::String = "NiwotRidge",
            g1::Bool = false
) where {FT<:AbstractFloat}
    @info tinfo("Plot main text example figure for $(site)...");

    # make sure that site is supported
    @assert site in ["NiwotRidge", "Ozark"];

    # use latex and serif
    if use_latex use_serif_tex(); end

    # replace this with artifact when all simulations are done
    file_loc = "/home/wyujie/RAID/Data/FLUXNET2015/simulation";
    #file_loc = artifact"2021_forest_tower_results";

    # read data as dataframe
    if g1
        df_bbm = read_csv("$(file_loc)/simulated_results_g_bbm_$(site).csv");
        df_med = read_csv("$(file_loc)/simulated_results_g_med_$(site).csv");
    else
        df_bbm = read_csv("$(file_loc)/simulated_results_bbm_$(site).csv");
        df_med = read_csv("$(file_loc)/simulated_results_med_$(site).csv");
    end
    df_osm = read_csv("$(file_loc)/simulated_results_osm_$(site).csv");
    df_bbm.YEAR = Int.(floor.(df_bbm.TIME ./ 1e8));
    df_med.YEAR = Int.(floor.(df_med.TIME ./ 1e8));
    df_osm.YEAR = Int.(floor.(df_osm.TIME ./ 1e8));
    df_bbm.DOY  = df_bbm.Day .+ (df_bbm.Hour .+ df_bbm.Minu ./ 60) ./ 24;
    df_med.DOY  = df_med.Day .+ (df_med.Hour .+ df_med.Minu ./ 60) ./ 24;
    df_osm.DOY  = df_osm.Day .+ (df_osm.Hour .+ df_osm.Minu ./ 60) ./ 24;

    # calculate the daily sum flow rates
    day_xmp = FT[];
    obs_CO₂ = FT[]; obs_H₂O = FT[];
    bbm_CO₂ = FT[]; bbm_H₂O = FT[];
    med_CO₂ = FT[]; med_H₂O = FT[];
    osm_CO₂ = FT[]; osm_H₂O = FT[];
    fac_CO₂ = 1e-6 * 0.5 * 3600 * 12 * -1;
    fac_H₂O = 0.5 * 3600 * M_H₂O();
    for iDay in 1:366
        mask = (abs.(df_bbm.Day .- iDay) .< 0.01) .* (df_osm.YEAR .== 2014);
        if sum(mask) > 0
            push!(day_xmp, iDay);
            push!(obs_CO₂, sum(df_bbm.ObsC[mask]) * fac_CO₂);
            push!(obs_H₂O, sum(df_bbm.ObsE[mask]) * fac_H₂O);
            push!(bbm_CO₂, sum(df_bbm.ModC[mask]) * fac_CO₂);
            push!(bbm_H₂O, sum(df_bbm.ModE[mask]) * fac_H₂O);
            push!(med_CO₂, sum(df_med.ModC[mask]) * fac_CO₂);
            push!(med_H₂O, sum(df_med.ModE[mask]) * fac_H₂O);
            push!(osm_CO₂, sum(df_osm.ModC[mask]) * fac_CO₂);
            push!(osm_H₂O, sum(df_osm.ModE[mask]) * fac_H₂O);
        end
    end

    # create canvas
    fig,axs = create_canvas("FT-example-$(site)"; nrow=4, figsize=(6,8.5));
    ax1,ax2,ax3,ax4 = axs;

    ax1.plot(day_xmp, obs_CO₂, "gray", label="obs");
    ax1.plot(day_xmp, bbm_CO₂, "tomato", label="BBM", alpha=0.6);
    ax1.plot(day_xmp, med_CO₂, "royalblue", label="MED", alpha=0.6);
    ax1.plot(day_xmp, osm_CO₂, "c", label="OSM", alpha=0.6);

    ax2.plot(day_xmp, obs_H₂O, "gray");
    ax2.plot(day_xmp, bbm_H₂O, "tomato", alpha=0.6);
    ax2.plot(day_xmp, med_H₂O, "royalblue", alpha=0.6);
    ax2.plot(day_xmp, osm_H₂O, "c", alpha=0.6);

    mask = (df_osm.TIME .>= 201409130000) .* (df_osm.TIME .<= 201409142359);

    ax3.plot(df_osm.DOY[mask], -df_osm.ObsC[mask], "gray");
    ax3.plot(df_osm.DOY[mask], -df_bbm.ModC[mask], "tomato");
    ax3.plot(df_osm.DOY[mask], -df_med.ModC[mask], "royalblue");
    ax3.plot(df_osm.DOY[mask], -df_osm.ModC[mask], "c");

    ax4.plot(df_osm.DOY[mask], df_osm.ObsE[mask] * 1000, "gray");
    ax4.plot(df_osm.DOY[mask], df_bbm.ModE[mask] * 1000, "tomato");
    ax4.plot(df_osm.DOY[mask], df_med.ModE[mask] * 1000, "royalblue");
    ax4.plot(df_osm.DOY[mask], df_osm.ModE[mask] * 1000, "c");

    ax1.legend(loc="lower right", ncol=4);
    set_xlims!(axs[3:4], [256,258]);
    set_xticks!(axs[3:4], [256,257,258]);
    set_xlabels!(axs, ["","","","Day of year 2014"]);
    set_ylabels!(axs, ["NEE (g C m\$^{-2}\$ day\$^{-1}\$)",
                       "ET (kg H\$_2\$O m\$^{-2}\$ day\$^{-1}\$)",
                       LS_NEE_unit, LS_ET_unit_m],
                 fontsize=12);
    set_titles!(axs; loc="left");

    save_canvas!(fig, "figures/2021_forest_tower/$(nTH)_example_$(site).pdf",
                 saving);
    if saving close(fig) end;

    return nothing
end








###############################################################################
#
# Figs. 9 and 10: Fitted variables for ForestTower2021 project
#
###############################################################################
"""
    plot_FT2021_fitting!(
                proj::ForestTower2021{FT},
                nTH::Int;
                saving::Bool = false,
                use_latex::Bool = true,
                site::String = "NiwotRidge"
    ) where {FT<:AbstractFloat}

Plot model fitting results for project ForestTower2021 for a site, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `nTH` Figure number
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
- `site` Site of flux tower, must be `NiwotRidge` or `Ozark`
"""
function plot_FT2021_fitting!(
            proj::ForestTower2021{FT},
            nTH::Int;
            saving::Bool = false,
            use_latex::Bool = true,
            site::String = "NiwotRidge"
) where {FT<:AbstractFloat}
    @info tinfo("Plot main text fitted variable figure for $(site)...");

    # make sure that site is supported
    @assert site in ["NiwotRidge", "Ozark"];

    # use latex and serif
    if use_latex use_serif_tex(); end

    # replace this with artifact when all simulations are done
    file_loc = "/home/wyujie/RAID/Data/FLUXNET2015/simulation";
    #file_loc = artifact"2021_forest_tower_results";

    # read data as dataframe
    df_bbm = read_csv("$(file_loc)/simulated_results_bbm_$(site).csv");
    df_med = read_csv("$(file_loc)/simulated_results_med_$(site).csv");
    df_osm = read_csv("$(file_loc)/simulated_results_osm_$(site).csv");
    df_bbm.YEAR = Int.(floor.(df_bbm.TIME ./ 1e8));
    df_med.YEAR = Int.(floor.(df_med.TIME ./ 1e8));
    df_osm.YEAR = Int.(floor.(df_osm.TIME ./ 1e8));

    # calculate the MAPEs
    mapec_bbm = []; mapec_med = []; mapec_osm = [];
    mapee_bbm = []; mapee_med = []; mapee_osm = [];
    mape_bbm  = []; mape_med  = []; mape_osm  = [];
    for year in 2000:2020
        _m = (df_bbm.YEAR .== year);
        if sum(_m) > 0
            stdc  = nanstd(df_bbm.ObsC[_m]);
            stde  = nanstd(df_bbm.ObsE[_m]);
            bbm_c = nanmean(abs.(df_bbm.ModC[_m] .- df_bbm.ObsC[_m])) / stdc;
            bbm_e = nanmean(abs.(df_bbm.ModE[_m] .- df_bbm.ObsE[_m])) / stde;
            med_c = nanmean(abs.(df_med.ModC[_m] .- df_med.ObsC[_m])) / stdc;
            med_e = nanmean(abs.(df_med.ModE[_m] .- df_med.ObsE[_m])) / stde;
            osm_c = nanmean(abs.(df_osm.ModC[_m] .- df_osm.ObsC[_m])) / stdc;
            osm_e = nanmean(abs.(df_osm.ModE[_m] .- df_osm.ObsE[_m])) / stde;
            push!(mapec_bbm, bbm_c);
            push!(mapee_bbm, bbm_e);
            push!(mapec_med, med_c);
            push!(mapee_med, med_e);
            push!(mapec_osm, osm_c);
            push!(mapee_osm, osm_e);
            push!(mape_bbm, bbm_c + bbm_e);
            push!(mape_med, med_c + med_e);
            push!(mape_osm, osm_c + osm_e);
        end
    end

    # plot fitting comparisons
    fit_bbm = read_csv("$(file_loc)/fitted_vrk_bbm_$(site).csv");
    fit_med = read_csv("$(file_loc)/fitted_vrk_med_$(site).csv");
    fit_osm = read_csv("$(file_loc)/fitted_vrk_osm_$(site).csv");
    fit_bbm.MAPEc = mapec_bbm;
    fit_med.MAPEc = mapec_med;
    fit_osm.MAPEc = mapec_osm;
    fit_bbm.MAPEe = mapee_bbm;
    fit_med.MAPEe = mapee_med;
    fit_osm.MAPEe = mapee_osm;
    fit_bbm.MAPE  = mape_bbm ;
    fit_med.MAPE  = mape_med ;
    fit_osm.MAPE  = mape_osm ;

    fig,axs = create_canvas("FT-fitting-$(site)"; nrow=2, figsize=(5,4.5));
    ax1,ax2 = axs;
    tx1     = ax1.twinx();

    ms = ["BBM", "MED", "OSM"];
    cs = ["tomato", "royalblue", "c"];
    ds = [fit_bbm, fit_med, fit_osm];

    # 300X for NiwotRidge and 100X for Ozark
    fk = site=="NiwotRidge" ? 400 : 100;

    for i in 1:3
        x1 = [0,4,8] .+ i;
        x2 = [0,4,8] .+ i;
        # filter out outliers
        if (site=="NiwotRidge") || (i==3)
            _m = ds[i].K .< 1;
        else
            _m = ds[i].K .> 0;
        end
        y1 = [nanmean(ds[i].Vcmax[ds[i].Vcmax .< 80]),
              nanmean(ds[i].Rbase[ds[i].Vcmax .< 80]),
              nanmean(ds[i].K[_m]) / 2 * fk];
        e1 = [nanstd(ds[i].Vcmax[ds[i].Vcmax .< 80]),
              nanstd(ds[i].Rbase[ds[i].Vcmax .< 80]),
              nanstd(ds[i].K[_m]) / 2 * fk];
        y2 = [nanmean(ds[i].MAPEc),
              nanmean(ds[i].MAPEe),
              nanmean(ds[i].MAPE)] .* 100;
        e2 = [nanstd(ds[i].MAPEc),
              nanstd(ds[i].MAPEe),
              nanstd(ds[i].MAPE)] .* 100;
        ax1.bar(x1, y1, width=0.9, yerr=e1, color=cs[i], label=ms[i]);
        ax2.bar(x2, y2, width=0.9, yerr=e2, color=cs[i], label=ms[i]);

        @show y1,e1;
        @show y2,e2;
    end

    t1 = [LS_Vcmax, LS_Rbase, LS_Kmax];
    t2 = ["NEE", "ET", "NEE\\&ET"];
    set_xticks!(axs, [[2,6,10],[2,6,10]]);
    set_ylabels!([ax1,ax2,tx1],
                 ["$(LS_Vcmax) or $(LS_Rbase)\n"*latex_unit("A"), "MASE (\\%)",
                  LS_Kmax*"\n(mol s\$^{-1}\$ MPa\$^{-1}\$)"]);
    if site == "NiwotRidge"
        set_ylims!([ax1,tx1], [[0,44], [0,0.11]]);
    else
        set_ylims!([ax1,tx1], [[0,60], [0,0.60]]);
        ax1.text(9.5, 60.2, "$(LS_Kmax) > 3 for BBM and MED",
                 ha="center", va="bottom");
    end
    ax1.set_xticklabels(t1, fontsize=16);
    ax2.set_xticklabels(t2, fontsize=16);
    ax1.legend(loc="upper center");
    set_titles!(axs; loc="left");

    save_canvas!(fig, "figures/2021_forest_tower/$(nTH)_fittings_$(site).pdf",
                 saving);
    if saving close(fig) end;

    return nothing
end








###############################################################################
#
# Figs. 11-14: Hex bin comparisons for ForestTower2021 project
#
###############################################################################
"""
    plot_FT2021_comparison!(
                proj::ForestTower2021{FT},
                nTH::Int;
                saving::Bool = false,
                use_latex::Bool = true,
                site::String = "NiwotRidge",
                g1::Bool = false
    ) where {FT<:AbstractFloat}

Plot model comparison hexbins for project ForestTower2021 for a site, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `nTH` Figure number
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
- `site` Site of flux tower, must be `NiwotRidge` or `Ozark`
- `g1` If `true`, empirical g1 was fitted. Default is `false`
"""
function plot_FT2021_comparison!(
            proj::ForestTower2021{FT},
            nTH::Int;
            saving::Bool = false,
            use_latex::Bool = true,
            site::String = "NiwotRidge",
            g1::Bool = false
) where {FT<:AbstractFloat}
    @info tinfo("Plot main text model comparison hexbin for $(site)...");

    # make sure that site is supported
    @assert site in ["NiwotRidge", "Ozark"];

    # use latex and serif
    if use_latex use_serif_tex(); end

    # replace this with artifact when all simulations are done
    file_loc = "/home/wyujie/RAID/Data/FLUXNET2015/simulation";
    #file_loc = artifact"2021_forest_tower_results";

    # read data as dataframe
    if g1
        df_bbm = read_csv("$(file_loc)/simulated_results_g_bbm_$(site).csv");
        df_med = read_csv("$(file_loc)/simulated_results_g_med_$(site).csv");
    else
        df_bbm = read_csv("$(file_loc)/simulated_results_bbm_$(site).csv");
        df_med = read_csv("$(file_loc)/simulated_results_med_$(site).csv");
    end
    df_osm = read_csv("$(file_loc)/simulated_results_osm_$(site).csv");

    # create canvas
    fig,axs = create_canvas("FT-density-$(site)"; figsize=(7,5.39),
                            ncol=3, nrow=2);
    ax1,ax2,ax3,ax4,ax5,ax6 = axs;

    # plot the hex bins
    if site=="NiwotRidge"
        c_lims,w_lims = [-13,22],[-3e-3,15e-3];
    else # elseif site=="Ozark"
        c_lims,w_lims = [-42,46],[-10e-3,21e-3];
    end
    plot_hexbin!(ax1, df_bbm.ObsC, df_bbm.ModC, c_lims, c_lims; logbins=true);
    plot_hexbin!(ax2, df_med.ObsC, df_med.ModC, c_lims, c_lims; logbins=true);
    plot_hexbin!(ax3, df_osm.ObsC, df_osm.ModC, c_lims, c_lims; logbins=true);
    plot_hexbin!(ax4, df_bbm.ObsE, df_bbm.ModE, w_lims, w_lims; logbins=true);
    plot_hexbin!(ax5, df_med.ObsE, df_med.ModE, w_lims, w_lims; logbins=true);
    plot_hexbin!(ax6, df_osm.ObsE, df_osm.ModE, w_lims, w_lims; logbins=true);

    # plot linear fitting curve and 1:1
    plot_line_regress!(ax1, df_bbm.ObsC, df_bbm.ModC);
    plot_line_regress!(ax2, df_med.ObsC, df_med.ModC);
    plot_line_regress!(ax3, df_osm.ObsC, df_osm.ModC);
    plot_line_regress!(ax4, df_bbm.ObsE, df_bbm.ModE);
    plot_line_regress!(ax5, df_med.ObsE, df_med.ModE);
    plot_line_regress!(ax6, df_osm.ObsE, df_osm.ModE);
    for ax in [ax1,ax2,ax3] ax.plot(c_lims, c_lims, "k:") end;
    for ax in [ax4,ax5,ax6] ax.plot(w_lims, w_lims, "k:") end;

    # display the fitting results table
    lrc_bbm = line_regress(df_bbm.ObsC, df_bbm.ModC);
    lrc_med = line_regress(df_med.ObsC, df_med.ModC);
    lrc_osm = line_regress(df_osm.ObsC, df_osm.ModC);
    lre_bbm = line_regress(df_bbm.ObsE, df_bbm.ModE);
    lre_med = line_regress(df_med.ObsE, df_med.ModE);
    lre_osm = line_regress(df_osm.ObsE, df_osm.ModE);
    lrcs = [lrc_bbm, lrc_med, lrc_osm];
    lres = [lre_bbm, lre_med, lre_osm];
    p_c_bbm = test_slope(df_bbm.ObsC, df_bbm.ModC; slope=1);
    p_c_med = test_slope(df_med.ObsC, df_med.ModC; slope=1);
    p_c_osm = test_slope(df_osm.ObsC, df_osm.ModC; slope=1);
    p_e_bbm = test_slope(df_bbm.ObsE, df_bbm.ModE; slope=1);
    p_e_med = test_slope(df_med.ObsE, df_med.ModE; slope=1);
    p_e_osm = test_slope(df_osm.ObsE, df_osm.ModE; slope=1);

    df_lrs = DataFrame();
    df_lrs.r2C    = [_lr.r2    for _lr in lrcs];
    df_lrs.interC = [_lr.inter for _lr in lrcs];
    df_lrs.slopeC = [_lr.slope for _lr in lrcs];
    df_lrs.pC     = [p_c_bbm, p_c_med, p_c_osm];
    df_lrs.r2E    = [_lr.r2    for _lr in lres];
    df_lrs.interE = [_lr.inter for _lr in lres];
    df_lrs.slopeE = [_lr.slope for _lr in lres];
    df_lrs.pE     = [p_e_bbm, p_e_med, p_e_osm];
    @show df_lrs;

    # set titles and labels
    set_yticklabels!([ax2,ax3,ax5,ax6], [[],[],[],[]]);
    _tlabs = ["BBM", "MED", "OSM", "BBM", "MED", "OSM"];
    _xlabs = ["", "obs \$-\$$(LS_NEE_unit)", "", "", "obs $(LS_ET_unit)", ""];
    _ylabs = ["mod \$-\$NEE\n" * latex_unit("A"), "", "",
              "mod ET\n" * latex_unit("E"), "", ""];
    set_titles!(axs; labels=_tlabs, loc="left");
    set_xylabels!(axs, _xlabs, _ylabs);
    set_xylims!(axs[1:3], c_lims, c_lims);
    set_xylims!(axs[4:6], w_lims, w_lims);
    set_titles!(axs; labels=_tlabs, loc="left");

    save_canvas!(fig, "figures/2021_forest_tower/$(nTH)_density_$(site).pdf",
                 saving);
    if saving close(fig) end;

    return nothing
end








###############################################################################
#
# Fig. 15: SIF comparisons for ForestTower2021 project
#
###############################################################################
"""
    plot_FT2021_SIF!(
                proj::ForestTower2021{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot SIF comparisons for project ForestTower2021, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_FT2021_SIF!(
            proj::ForestTower2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    @info tinfo("Plot main text SIF comparisons...");

    # use latex and serif
    if use_latex use_serif_tex(); end

    # replace this with artifact when all simulations are done
    file_loc = "/home/wyujie/RAID/Data/FLUXNET2015/simulation";
    #file_loc = artifact"2021_forest_tower_results";

    # read gridded TROPOMI SIF

    # simulated SIF
    df_2018n = read_csv("$(file_loc)/sif_NiwotRidge_2018.csv");
    df_2019n = read_csv("$(file_loc)/sif_NiwotRidge_2019.csv");
    df_nr1   = vcat([df_2018n, df_2019n]...);
    df_2018o = read_csv("$(file_loc)/sif_Ozark_2018.csv");
    df_2019o = read_csv("$(file_loc)/sif_Ozark_2019.csv");
    df_moz   = vcat([df_2018o, df_2019o]...);
    lf_2018n = read_csv("$(file_loc)/sif_NiwotRidge_2018_LAI.csv");
    lf_2019n = read_csv("$(file_loc)/sif_NiwotRidge_2019_LAI.csv");
    lf_nr1   = vcat([lf_2018n, lf_2019n]...);
    lf_2018o = read_csv("$(file_loc)/sif_Ozark_2018_LAI.csv");
    lf_2019o = read_csv("$(file_loc)/sif_Ozark_2019_LAI.csv");
    lf_moz   = vcat([lf_2018o, lf_2019o]...);

    # plot figures
    fig,axs = create_canvas("FT-SIF"; ncol=2, figsize=(6,3.5));
    ax1,ax2 = axs;

    # mask the data by Coverage
    mask_nr1 = (df_nr1.Coverage .>= 50) .* (df_nr1.obsSIF .< 2);
    mask_moz = df_moz.Coverage .>= 80;

    ax1.plot(df_nr1.obsSIF[mask_nr1], df_nr1.SIF740[mask_nr1], "ko";
             label="Site LAI");
    ax1.plot(lf_nr1.obsSIF[mask_nr1], lf_nr1.SIF740[mask_nr1], "c+";
             label="MODIS LAI");
    ax2.plot(df_moz.obsSIF[mask_moz], df_moz.SIF740[mask_moz], "ko";
             label="Site LAI");
    ax2.plot(lf_moz.obsSIF[mask_moz], lf_moz.SIF740[mask_moz], "c+";
             label="MODIS LAI");
    plot_line_regress!(ax1, df_nr1.obsSIF[mask_nr1], df_nr1.SIF740[mask_nr1];
                       interval=true, color="k");
    plot_line_regress!(ax1, lf_nr1.obsSIF[mask_nr1], lf_nr1.SIF740[mask_nr1];
                       interval=true, color="c");
    plot_line_regress!(ax2, df_moz.obsSIF[mask_moz], df_moz.SIF740[mask_moz];
                       interval=true, color="k");
    plot_line_regress!(ax2, lf_moz.obsSIF[mask_moz], lf_moz.SIF740[mask_moz];
                       interval=true, color="c");

    ax1.plot([-1.0,2.0], [-1.0,2.0], "k:");
    ax2.plot([-0.5,3.5], [-0.5,3.5], "k:");

    rmse_11 = round(rmse(df_nr1.obsSIF[mask_nr1], df_nr1.SIF740[mask_nr1]);
                    digits=3);
    rmse_12 = round(rmse(lf_nr1.obsSIF[mask_nr1], lf_nr1.SIF740[mask_nr1]);
                    digits=3);
    rmse_21 = round(rmse(df_moz.obsSIF[mask_moz], df_moz.SIF740[mask_moz]);
                    digits=3);
    rmse_22 = round(rmse(lf_moz.obsSIF[mask_moz], lf_moz.SIF740[mask_moz]);
                    digits=3);
    ax1.text(2.0, -0.8, "RMSE = $(rmse_11)", ha="right", color="k",
             va="bottom", fontsize=12);
    ax1.text(2.0, -0.8, "RMSE = $(rmse_12)", ha="right", color="c",
             va="top", fontsize=12);
    ax2.text(3.5, -0.3, "RMSE = $(rmse_21)", ha="right", color="k",
             va="bottom", fontsize=12);
    ax2.text(3.5, -0.3, "RMSE = $(rmse_22)", ha="right", color="c",
             va="top", fontsize=12);

    lr_n = line_regress(df_nr1.obsSIF[mask_nr1], df_nr1.SIF740[mask_nr1]);
    @show lr_n;
    lr_n = line_regress(lf_nr1.obsSIF[mask_nr1], lf_nr1.SIF740[mask_nr1]);
    @show lr_n;
    lr_o = line_regress(df_moz.obsSIF[mask_moz], df_moz.SIF740[mask_moz]);
    @show lr_o;
    lr_o = line_regress(lf_moz.obsSIF[mask_moz], lf_moz.SIF740[mask_moz]);
    @show lr_o;

    ax1.set_yticks([-1,0,1,2])
    ax2.legend(loc="upper left");
    set_xlabels!(axs, "obs SIF\$_{740}\$\n$(LS_unit_SIF)");
    set_ylabels!([ax1], "OSM SIF\$_{740}\$\n$(LS_unit_SIF)");
    set_titles!(axs; loc="left");

    save_canvas!(fig, "figures/2021_forest_tower/15_sif_series.pdf", saving);
    if saving close(fig) end;

    return nothing;
end








###############################################################################
#
# Plot all main text figures for ForestTower2021 project
#
###############################################################################
"""
    plot_FT2021!(
                proj::ForestTower2021{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot figures for project ForestTower2021, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_FT2021!(
            proj::ForestTower2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    plot_FT2021_ecrit!(proj; saving=saving, use_latex=use_latex);
    plot_FT2021_spectrum!(proj; saving=saving, use_latex=use_latex);
    plot_FT2021_clumping!(proj; saving=saving, use_latex=use_latex);
    plot_FT2021_responses!(proj; saving=saving, use_latex=use_latex);
    plot_FT2021_example!(proj, 7; saving=saving, use_latex=use_latex,
                         site="NiwotRidge");
    plot_FT2021_example!(proj, 8; saving=saving, use_latex=use_latex,
                         site="Ozark");
    plot_FT2021_fitting!(proj, 9; saving=saving, use_latex=use_latex,
                         site="NiwotRidge");
    plot_FT2021_fitting!(proj, 10; saving=saving, use_latex=use_latex,
                         site="Ozark");
    plot_FT2021_comparison!(proj, 11; saving=saving, use_latex=use_latex,
                            site="NiwotRidge");
    plot_FT2021_comparison!(proj, 12; saving=saving, use_latex=use_latex,
                            site="Ozark");
    plot_FT2021_comparison!(proj, 13; saving=saving, use_latex=use_latex,
                            site="NiwotRidge", g1=true);
    plot_FT2021_comparison!(proj, 14; saving=saving, use_latex=use_latex,
                            site="Ozark", g1=true);
    plot_FT2021_SIF!(proj; saving=saving, use_latex=use_latex);

    return nothing;
end








###############################################################################
#
# Plot all SI figures
#
###############################################################################
"""
    plot_FT2021_SI!(
                proj::ForestTower2021{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot figures for project ForestTower2021, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_FT2021_SI!(
            proj::ForestTower2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    # plot comparison for each year
    plot_FT2021_SI_example!(proj; saving=saving, use_latex=use_latex,
                            site="NiwotRidge");
    plot_FT2021_SI_example!(proj; saving=saving, use_latex=use_latex,
                            site="Ozark");

    return nothing;
end








###############################################################################
#
# Plot testing figures, not for main text purpose
#
###############################################################################
"""
    plot_FT2021_SI_example!(
                proj::ForestTower2021{FT};
                saving::Bool = false,
                use_latex::Bool = true,
                site::String = "NiwotRidge
    ) where {FT<:AbstractFloat}

Plot SI figures for project ForestTower2021 for each year, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
- `site` Site of flux tower, must be `NiwotRidge` or `Ozark`
"""
function plot_FT2021_SI_example!(
            proj::ForestTower2021{FT};
            saving::Bool = false,
            use_latex::Bool = true,
            site::String = "NiwotRidge"
) where {FT<:AbstractFloat}
    @info tinfo("Plot SI figures for each year at $(site)...");

    # make sure that site is supported
    @assert site in ["NiwotRidge", "Ozark"];

    # use latex and serif
    if use_latex use_serif_tex(); end

    # replace this with artifact when all simulations are done
    file_loc = "/home/wyujie/RAID/Data/FLUXNET2015/simulation";
    #file_loc = artifact"2021_forest_tower_results";

    # read data as dataframe
    df_bbm = read_csv("$(file_loc)/simulated_results_bbm_$(site).csv");
    df_med = read_csv("$(file_loc)/simulated_results_med_$(site).csv");
    df_osm = read_csv("$(file_loc)/simulated_results_osm_$(site).csv");

    df_bbm.YEAR = Int.(floor.(df_bbm.TIME ./ 1e8));
    df_med.YEAR = Int.(floor.(df_med.TIME ./ 1e8));
    df_osm.YEAR = Int.(floor.(df_osm.TIME ./ 1e8));
    df_bbm.DOY  = df_bbm.Day .+ (df_bbm.Hour .+ df_bbm.Minu ./ 60) ./ 24;
    df_med.DOY  = df_med.Day .+ (df_med.Hour .+ df_med.Minu ./ 60) ./ 24;
    df_osm.DOY  = df_osm.Day .+ (df_osm.Hour .+ df_osm.Minu ./ 60) ./ 24;

    # create canvas
    for year in 2006:2019
        mask = (df_osm.YEAR .== year);

        # if there are enough data
        if sum(mask) > 2400
            fig,axs = create_canvas("FT-$(site)-$(year)"; nrow=2,
                                    figsize=(24,5));
            ax1,ax2 = axs;

            ax1.plot(df_osm.DOY[mask], df_osm.ObsC[mask], "gray", linewidth=1);
            ax1.plot(df_bbm.DOY[mask], df_bbm.ModC[mask], "tomato",
                     alpha=0.6);
            ax1.plot(df_med.DOY[mask], df_med.ModC[mask], "royalblue",
                     alpha=0.6);
            ax1.plot(df_osm.DOY[mask], df_osm.ModC[mask], "c",
                     alpha=0.6);

            ax2.plot(df_osm.DOY[mask], df_osm.ObsE[mask], "gray", linewidth=1);
            ax2.plot(df_bbm.DOY[mask], df_bbm.ModE[mask], "tomato",
                     alpha=0.6);
            ax2.plot(df_med.DOY[mask], df_med.ModE[mask], "royalblue",
                     alpha=0.6);
            ax2.plot(df_osm.DOY[mask], df_osm.ModE[mask], "c",
                     alpha=0.6);

            set_titles!(axs; loc="left");
            set_xlabels!(axs, ["","Day of year $(year)"]);
            set_ylabels!(axs, ["\$-\$"*LS_NEE_unit, LS_ET_unit], fontsize=12);

            fn = "figures/2021_forest_tower/per_year/SI-$(site)-$(year).pdf";
            save_canvas!(fig, fn, saving);
            if saving close(fig) end;
        end
    end

    return nothing;
end




"""
    plot_FT2021_patm!(
                proj::ForestTower2021{FT};
                saving::Bool = false,
                use_latex::Bool = true,
                site::String = "NiwotRidge
    ) where {FT<:AbstractFloat}

Plot test figures for project ForestTower2021 at different atmospheric pressure
    settings, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_FT2021_patm!(
            proj::ForestTower2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    @info tinfo("Plot test figures for P_atm...");

    # use latex and serif
    if use_latex use_serif_tex(); end

    # replace this with artifact when all simulations are done
    file_loc = "/home/wyujie/RAID/Data/FLUXNET2015/simulation";
    #file_loc = artifact"2021_forest_tower_results";

    # read data
    df_patm = read_csv("$(file_loc)/df_flux_tower.csv");
    df_seal = read_csv("$(file_loc)/df_sea_level.csv" );

    # plot the results
    fig,axs = create_canvas("FT-patm"; nrow=2, figsize=(24,5));
    ax1,ax2 = axs;

    ax1.plot(df_patm.DOY, df_patm.ObsC, "gray", label="obs", linewidth=1);
    ax1.plot(df_patm.DOY, df_patm.ModC, "c", alpha=0.6, label="atm",
             linewidth=1);
    ax1.plot(df_patm.DOY, df_seal.ModC, "tomato", alpha=0.6, label="sea",
             linewidth=1);
    ax1.legend(loc="lower right", ncol=3);

    ax2.plot(df_patm.DOY, df_patm.ObsE, "gray", label="obs", linewidth=1);
    ax2.plot(df_patm.DOY, df_patm.ModE, "c", alpha=0.6, label="atm",
             linewidth=1);
    ax2.plot(df_patm.DOY, df_seal.ModE, "tomato", alpha=0.6, label="sea",
             linewidth=1);

    set_titles!(axs; loc="left");
    set_xlabels!(axs, ["","Day of year 2019"]);
    set_ylabels!(axs, [LS_NEE_unit, LS_ET_unit], fontsize=12);

    save_canvas!(fig, "figures/2021_forest_tower/SI-patm.pdf", saving);
    if saving close(fig) end;

    return nothing;
end
