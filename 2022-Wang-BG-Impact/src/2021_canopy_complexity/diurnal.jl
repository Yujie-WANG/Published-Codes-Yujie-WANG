###############################################################################
#
# Test how the methods differ in simulated diurnal cycle
#
###############################################################################
"""
This function simulates the diurnal cycles of the SPACMono node using different
    canopy layering modes.
"""
function diurnal_cycles! end




"""
This method simulates the diurnal cycles and then plot the results:
    diurnal_cycles!(
                proj::CanopyComplexity2021{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Simulate and plot the diurnal cycles of different canopy layering modes, given
- `proj` [`CanopyComplexity2021`](@ref) type project identifier
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
diurnal_cycles!(
            proj::CanopyComplexity2021{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat} = (
    # use latex and serif
    if use_latex use_serif_tex(); end;
    @info tinfo("Plot diurnal cycles...");

    # create nodes to work on
    _node_ijkx = create_spac(proj, "IJKX");
    _node_2kx  = create_spac(proj, "2KX");
    _node_kx   = create_spac(proj, "KX");
    _node_2x   = create_spac(proj, "2X");
    _node_1x   = create_spac(proj, "1X");
    _nodes = [_node_ijkx, _node_2kx, _node_kx, _node_2x, _node_1x];
    _soil_hs = VanGenuchten{FT}(stype = "Ozark", α = 1.368, n = 2.6257, Θs = 0.45, Θr = 0.067);
    for _node in _nodes
        for _root in _node.plant_hs.roots
            _root.sh = deepcopy(_soil_hs);
        end;
        _node.ga                = 413.223;
        _node.la                = 1735.537;
        _node.latitude          = 38.74;
        _node.longitude         = -92.2;
        _node.elevation         = 219.4;
        _node.canopy_rt.Ω       = 0.69;
        _node.canopy_rt.clump_a = 0.69;
        for _iPS in _node.plant_ps
            _iPS.g_min   = 0.001;
            _iPS.g_min25 = 0.001;
            _iPS.g_max   = 0.095;
            _iPS.g_max25 = 0.095;
        end;
        update_LAI!(_node, FT(4.2));
        update_VJRWW!(_node, FT(75));
        update_Weibull!(_node, FT(5.703), FT(0.953));
        update_Kmax!(_node, FT(1.55));
        update_Cab!(_node, FT(57.23));
        mod_spac!(_node_ijkx, "C3");
        mod_spac!(_node_2kx, "C3");
        mod_spac!(_node_kx, "C3");
    end;

    # read Ozarks data and create new columns
    _rawdf = query_data(proj, 2019, "Ozark");
    _days  = [177,180];
    _newdf = _rawdf[_days[1] .<= _rawdf.Day .< _days[2],:];
    _time  = _newdf.Day .+ (_newdf.Hour .+ (_newdf.Minu .+ 15) ./  60) ./ 24;
    _newdf.C_IJKX = [NaN for _i in eachindex(_newdf.Day)];
    _newdf.C_2KX  = [NaN for _i in eachindex(_newdf.Day)];
    _newdf.C_KX   = [NaN for _i in eachindex(_newdf.Day)];
    _newdf.C_2X   = [NaN for _i in eachindex(_newdf.Day)];
    _newdf.C_1X   = [NaN for _i in eachindex(_newdf.Day)];
    _newdf.E_IJKX = [NaN for _i in eachindex(_newdf.Day)];
    _newdf.E_2KX  = [NaN for _i in eachindex(_newdf.Day)];
    _newdf.E_KX   = [NaN for _i in eachindex(_newdf.Day)];
    _newdf.E_2X   = [NaN for _i in eachindex(_newdf.Day)];
    _newdf.E_1X   = [NaN for _i in eachindex(_newdf.Day)];

    # run diurnal cycles
    diurnal_cycles!(proj, _nodes, _newdf);

    # change units
    _newdf.C_IJKX .*= -1;
    _newdf.C_2KX  .*= -1;
    _newdf.C_KX   .*= -1;
    _newdf.C_2X   .*= -1;
    _newdf.C_1X   .*= -1;
    _newdf.ObsC   .*= -1;
    _newdf.E_IJKX .*= 1000;
    _newdf.E_2KX  .*= 1000;
    _newdf.E_KX   .*= 1000;
    _newdf.E_2X   .*= 1000;
    _newdf.E_1X   .*= 1000;
    _newdf.ObsE   .*= 1000;

    # plot the simulations out
    _legends = ["IJKX", "2KX", "KX", "2X", "1X"];
    _datas_c = [_newdf.C_IJKX, _newdf.C_2KX, _newdf.C_KX, _newdf.C_2X, _newdf.C_1X];
    _datas_e = [_newdf.E_IJKX, _newdf.E_2KX, _newdf.E_KX, _newdf.E_2X, _newdf.E_1X];
    (_fig,(_ax1,_ax2)) = create_canvas("CC-diurnal"; nrow=2, figsize=(6,6));
    for _i in eachindex(_legends)
        _ax1.plot(_time, _datas_c[_i], label=_legends[_i], alpha=0.7);
        _ax2.plot(_time, _datas_e[_i], label=_legends[_i], alpha=0.7);
    end;
    _ax1.plot(_time, _newdf.ObsC, "k:", label="obs", alpha=0.7);
    _ax2.plot(_time, _newdf.ObsE, "k:", alpha=0.7);

    _ax1.legend(loc="upper center", ncol=6);
    set_xticks!([_ax1,_ax2], [177,178,179,180]);
    set_titles!([_ax1,_ax2]; loc="left");
    set_xylabels!([_ax1,_ax2], "Day of Year 2019", [LS_NEE_unit,LS_ET_unit_m]);
    save_canvas!(_fig, "figures/2021_canopy_complexity/7-diurnal.pdf", saving);
    if saving close(_fig) end;

    # plot the 1:1 comparisons
    (_fig,(_ax1,_ax2)) = create_canvas("CC-diurnal_comparison"; ncol=2);
    _mean_c = [nanmean(_dat .- _newdf.ObsC) for _dat in _datas_c];
    _mean_e = [nanmean(_dat .- _newdf.ObsE) for _dat in _datas_e];
    _stdv_c = [nanstd(_dat .- _newdf.ObsC) for _dat in _datas_c];
    _stdv_e = [nanstd(_dat .- _newdf.ObsE) for _dat in _datas_e];
    _ax1.bar([1,2,3,4,5], _mean_c, yerr=_stdv_c, color=COLORS[1:5], alpha=0.7);
    _ax2.bar([1,2,3,4,5], _mean_e, yerr=_stdv_e, color=COLORS[1:5], alpha=0.7);
    set_xticks!([_ax1,_ax2], [1,2,3,4,5]);
    _ax1.set_xticklabels(_legends, rotation=30);
    _ax2.set_xticklabels(_legends, rotation=30);
    set_titles!([_ax1,_ax2]; loc="left");
    set_ylabels!([_ax1,_ax2], [LS_NEE_unit,LS_ET_unit_m]);
    save_canvas!(_fig, "figures/2021_canopy_complexity/8-compare.pdf", saving);
    if saving close(_fig) end;

    return nothing
)




"""
This method simulates the diurnal cycles:
    diurnal_cycles!(
                proj::CanopyComplexity2021{FT},
                nodes::Vector{SPACMono{FT}},
                df::DataFrame,
                rbase::Q10TD{FT} = Q10TD{FT}(0, 298.15, 1.7)
    ) where {FT<:AbstractFloat}

Simulate and plot the diurnal cycles of different canopy layering modes, given
- `proj` [`CanopyComplexity2021`](@ref) type project identifier
- `nodes` Vector of 5 nodes with different canopy layering modes
- `df` Given data frame that stores the weather information
- `rbase` Base respiration rate for wood and soil
"""
diurnal_cycles!(
            proj::CanopyComplexity2021{FT},
            nodes::Vector{SPACMono{FT}},
            df::DataFrame,
            rbase::Q10TD{FT} = Q10TD{FT}(2, 298.15, 1.7)
) where {FT<:AbstractFloat} = (
    # use IJKX mode as default node
    _node_ijkx  = nodes[1];
    _in_rad_bak = deepcopy(_node_ijkx.in_rad);

    # 0.1 unpack data

    # iterate through the weather data
    @showprogress for _i in eachindex(df.Day)
        # update soil water matrices for all nodes
        _w_soil = FT(df.SWC[_i]) / 100;
        _t_soil = FT(df.T_SOIL[_i] + 273.15);
        _ψ_soil = soil_p_25_swc(_node_ijkx.plant_hs.roots[1].sh, _w_soil);
        for _node in nodes
            for _root in _node.plant_hs.roots
                _root.p_ups = _ψ_soil;
            end;
        end;

        # update PAR related information
        _zenith = zenith_angle(_node_ijkx.latitude, FT(df.Day[_i]), FT(df.Hour[_i]), FT(df.Minu[_i]));
        _zenith = min(88, _zenith);
        for _node in nodes
            mod_spac!(_node, _zenith);
            (_node).angles.sza = _zenith;
        end;

        # use
        _in_SW = (numerical∫(_in_rad_bak.E_direct, _node_ijkx.wl_set.dWL) + numerical∫(_in_rad_bak.E_diffuse, _node_ijkx.wl_set.dWL)) / 1000;
        _ratio = df.IN_RAD[_i] ./ _in_SW;
        # @show _in_rad_bak.E_direct;
        # @show _in_rad_bak.E_diffuse;
        # @show _in_SW, df.IN_RAD[_i],_ratio;
        # sleep(5);
        for _node in nodes
            _node.in_rad.E_direct  .= _in_rad_bak.E_direct  .* _ratio;
            _node.in_rad.E_diffuse .= _in_rad_bak.E_diffuse .* _ratio;
        end;

        # create an environment struct to use in all modes
        _envir = deepcopy(_node_ijkx.envirs[1]);
        _envir.t_air = df.T_AIR[_i] + 273.15;
        _envir.p_atm = df.P_ATM[_i] * 1000;
        _envir.p_a   = _envir.p_atm * 4e-4;
        _envir.p_O₂  = _envir.p_atm * 0.209;
        _envir.p_sat = saturation_vapor_pressure(_envir.t_air);
        _envir.vpd   = df.VPD[_i] * 100;
        _envir.p_H₂O = _envir.p_sat - _envir.vpd;
        _envir.RH    = _envir.p_H₂O / _envir.p_sat;
        _envir.wind  = df.WIND[_i];

        # prescribe leaf temperature
        _tl = (df.LW_OUT[_i] / 0.97 / K_STEFAN() ) ^ 0.25;
        for _node in nodes
            for _iPS in _node.plant_ps
                _iPS.T = _tl;
            end;
        end;

        # simulate the fluxes at steady state
        _nodes = deepcopy(nodes);
        simulation!(proj, _nodes, _envir; τ=FT(1e-5));

        # for debugging
        #=
        @show df.Hour[_i];
        @show nodes[4].plant_ps[1].APAR;
        @show nodes[4].plant_ps[2].APAR;
        @show nodes[4].plant_ps[1].LAIx;
        @show nodes[4].plant_ps[2].LAIx;
        @show nodes[4].plant_ps[1].LA;
        @show nodes[4].plant_ps[2].LA;
        @show nodes[4].plant_ps[1].ps.Vcmax25;
        @show nodes[4].plant_ps[2].ps.Vcmax25;
        @show nodes[4].plant_ps[1].ps.Vcmax;
        @show nodes[4].plant_ps[2].ps.Vcmax;
        @show nodes[5].plant_ps[1].APAR;
        @show nodes[5].plant_ps[1].LAIx;
        @show nodes[5].plant_ps[1].LA;
        @show nodes[5].plant_ps[1].ps.Vcmax25;
        @show nodes[5].plant_ps[1].ps.Vcmax;
        sleep(5);
        =#

        # calculate respiration from other components
        _r = photo_TD_from_set(rbase, _t_soil);

        # save data in the df
        df.C_IJKX[_i] = _nodes[1].f_npp - _r;
        df.C_2KX[_i]  = _nodes[2].f_npp - _r;
        df.C_KX[_i]   = _nodes[3].f_npp - _r;
        df.C_2X[_i]   = _nodes[4].f_npp - _r;
        df.C_1X[_i]   = _nodes[5].f_npp - _r;
        df.E_IJKX[_i] = _nodes[1].f_H₂O;
        df.E_2KX[_i]  = _nodes[2].f_H₂O;
        df.E_KX[_i]   = _nodes[3].f_H₂O;
        df.E_2X[_i]   = _nodes[4].f_H₂O;
        df.E_1X[_i]   = _nodes[5].f_H₂O;

        # update the data frame
        #df.ModC[_i] = _f_CO₂ / ga - _r;
        #df.ModE[_i] = _f_H₂O / ga;
    end;

    return nothing
)
