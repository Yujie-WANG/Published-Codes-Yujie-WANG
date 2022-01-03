###############################################################################
#
# Plot the spectra comparison
#
###############################################################################
"""
    spectra_diff!(
                proj::CanopyComplexity2021{FT},
                fname::String,
                hour::Int;
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot CanopyComplexity2021 SIF spectra comparison, given
- `proj` [`CanopyComplexity2021`](@ref) type project identifier
- `fname` Figure name to save
- `hour` Which time to simulate
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function spectra_diff!(
            proj::CanopyComplexity2021{FT},
            fname::String,
            hour::Int;
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    # use latex and serif
    if use_latex use_serif_tex(); end

    # create nodes and envir to work on
    @info tinfo("Plot BRDF effects on ΔSIF...");
    _nodes, _envir = create_spac(proj, OSMWang{FT}(), true, hour);
    simulation!(proj, _nodes, _envir; τ=FT(4e-6));

    # matrices for SIF at different VZA and RAA
    _vzas     = collect(FT, 0:5:85);
    _raas     = collect(FT, 0:10:360);
    _mat_ijkx = zeros(FT, (length(_vzas), length(_raas)));
    _mat_2kx  = zeros(FT, (length(_vzas), length(_raas)));
    _mat_kx   = zeros(FT, (length(_vzas), length(_raas)));
    _mat_2x   = zeros(FT, (length(_vzas), length(_raas)));
    _mat_1x   = zeros(FT, (length(_vzas), length(_raas)));
    @showprogress for _i in eachindex(_raas), _j in eachindex(_vzas)
        _mat_ijkx[_j,_i],_mat_2kx[_j,_i],_mat_kx[_j,_i],_mat_2x[_j,_i],_mat_1x[_j,_i] = spectra_diff!(proj, _nodes, _vzas[_j], _raas[_i]);
    end;

    # plot the data
    # TODO add an option of polar projection in PlotPlants.jl
    _fig = figure("CC-spectra-$(hour)"; figsize=(7.5,6.5), dpi=100);
    _fig.clear();
    _ax1 = _fig.add_subplot(2, 2, 1, polar=true);
    _ax2 = _fig.add_subplot(2, 2, 2, polar=true);
    _ax3 = _fig.add_subplot(2, 2, 3, polar=true);
    _ax4 = _fig.add_subplot(2, 2, 4, polar=true);

    _cm1 = _ax1.contourf(deg2rad.(_raas), _vzas, _mat_2kx ./ _mat_ijkx, levels=collect(1:0.01:1.2), extend="both");
    _cm2 = _ax2.contourf(deg2rad.(_raas), _vzas, _mat_kx  ./ _mat_ijkx, levels=collect(1:0.01:1.2), extend="both");
    _cm3 = _ax3.contourf(deg2rad.(_raas), _vzas, _mat_2x  ./ _mat_ijkx, levels=collect(1:0.01:1.2), extend="both");
    _cm4 = _ax4.contourf(deg2rad.(_raas), _vzas, _mat_1x  ./ _mat_ijkx, levels=collect(0.78:0.01:1.42), cmap="terrain_r");

    _ax1.set_yticks([]);
    _ax2.set_yticks([]);
    _ax3.set_yticks([]);
    _ax4.set_yticks([]);

    _ax1.grid(linestyle=":");
    _ax2.grid(linestyle=":");
    _ax3.grid(linestyle=":");
    _ax4.grid(linestyle=":");

    _fig.colorbar(_cm1; pad=0.08, ax=_ax1, ticks=collect(1.0:0.02:1.2));
    _fig.colorbar(_cm2; pad=0.08, ax=_ax2, ticks=collect(1.0:0.02:1.2));
    _fig.colorbar(_cm3; pad=0.08, ax=_ax3, ticks=collect(1.0:0.02:1.2));
    _fig.colorbar(_cm4; pad=0.08, ax=_ax4);

    set_titles!([_ax1,_ax2,_ax3,_ax4]; loc="left", labels=["2KX/IJKX", "KX/IJKX", "2X/IJKX", "1X/IJKX"]);
    save_canvas!(_fig, "figures/2021_canopy_complexity/$(fname).pdf", saving);
    if saving close(_fig) end;

    return nothing
end




"""
    spectra_diff!(
                proj::CanopyComplexity2021{FT},
                nodes::Vector{SPACMono{FT}},
                vza::FT,
                raa::FT
    ) where {FT<:AbstractFloat}

Plot CanopyComplexity2021 SIF spectra comparison, given
- `proj` [`CanopyComplexity2021`](@ref) type project identifier
- `nodes` Vector of 5 nodes with different canopy layering modes
- `vza` Viewing zenith angle
- `raa` Relative azimuth angle
"""
function spectra_diff!(
            proj::CanopyComplexity2021{FT},
            nodes::Vector{SPACMono{FT}},
            vza::FT,
            raa::FT
) where {FT<:AbstractFloat}
    # update environmental coonditions
    _node_ijkx, _node_2kx, _node_kx, _node_2x, _node_1x = nodes;
    _node_ijkx.angles.vza = vza;
    _node_ijkx.angles.raa = raa;

    # update canopy RT for IJKX node only
    for _node in nodes[1:1]
        @unpack angles, can_opt, can_rad, canopy_rt, in_rad, leaves_rt, rt_con, soil_opt, wl_set = _node;
        fluspect!.(leaves_rt, [wl_set]);
        canopy_geometry!(canopy_rt, angles, can_opt, rt_con);
        canopy_matrices!(leaves_rt, can_opt);
        short_wave!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, rt_con);
        canopy_fluxes!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, leaves_rt, wl_set, rt_con);
    end

    # update the PAR per layer for IJKX, 2KX, and KX modes
    @unpack can_opt, can_rad, canopy_rt, f_SL = _node_ijkx;
    @unpack nAzi, nIncl = canopy_rt;
    _nSL = nAzi * nIncl;

    # update fluorescence quantum yield for IJKX mode
    __node = deepcopy(_node_ijkx);
    for _i_can in 1:__node.n_canopy
        _iRT = __node.n_canopy + 1 - _i_can;
        _1PS = _node_ijkx.plant_ps[_i_can];
        __node.can_rad.ϕ_sun[:,:,_iRT] .= reshape(view((_1PS).φs,1:_nSL), nIncl, nAzi);
        __node.can_rad.ϕ_shade[_iRT] = (_1PS).φs[end];
    end
    SIF_fluxes!(__node.leaves_rt, __node.can_opt, __node.can_rad, __node.canopy_rt, __node.soil_opt, __node.wl_set, __node.rt_con, __node.rt_dim);
    _node_ijkx.can_rad = deepcopy(__node.can_rad);

    # update fluorescence quantum yield for 2KX mode
    __node = deepcopy(_node_ijkx);
    for _i_can in 1:__node.n_canopy
        _iRT = __node.n_canopy + 1 - _i_can;
        _2PS = _node_2kx.plant_ps[_i_can];
        __node.can_rad.ϕ_sun[:,:,_iRT] .= (_2PS).φs[1];
        __node.can_rad.ϕ_shade[_iRT] = (_2PS).φs[2];
    end
    SIF_fluxes!(__node.leaves_rt, __node.can_opt, __node.can_rad, __node.canopy_rt, __node.soil_opt, __node.wl_set, __node.rt_con, __node.rt_dim);
    _node_2kx.can_rad = deepcopy(__node.can_rad);

    # update fluorescence quantum yield for KX mode
    __node = deepcopy(_node_ijkx);
    for _i_can in 1:__node.n_canopy
        _iRT = __node.n_canopy + 1 - _i_can;
        _3PS = _node_kx.plant_ps[_i_can];
        __node.can_rad.ϕ_sun[:,:,_iRT] .= (_3PS).φs[1];
        __node.can_rad.ϕ_shade[_iRT] = (_3PS).φs[1];
    end
    SIF_fluxes!(__node.leaves_rt, __node.can_opt, __node.can_rad, __node.canopy_rt, __node.soil_opt, __node.wl_set, __node.rt_con, __node.rt_dim);
    _node_kx.can_rad = deepcopy(__node.can_rad);

    # update fluorescence quantum yield for 2X mode
    __node = deepcopy(_node_ijkx);
    for _i_can in 1:__node.n_canopy
        _iRT = __node.n_canopy + 1 - _i_can;
        _4PS = _node_2x.plant_ps;
        __node.can_rad.ϕ_sun[:,:,_iRT] .= (_4PS[1]).φs[1];
        __node.can_rad.ϕ_shade[_iRT] = (_4PS[2]).φs[1];
    end
    SIF_fluxes!(__node.leaves_rt, __node.can_opt, __node.can_rad, __node.canopy_rt, __node.soil_opt, __node.wl_set, __node.rt_con, __node.rt_dim);
    _node_2x.can_rad = deepcopy(__node.can_rad);

    # update fluorescence quantum yield for 1X mode
    __node = deepcopy(_node_ijkx);
    for _i_can in 1:__node.n_canopy
        _iRT = __node.n_canopy + 1 - _i_can;
        _5PS = _node_1x.plant_ps[1];
        __node.can_rad.ϕ_sun[:,:,_iRT] .= (_5PS).φs[1];
        __node.can_rad.ϕ_shade[_iRT] = (_5PS).φs[1];
    end
    SIF_fluxes!(__node.leaves_rt, __node.can_opt, __node.can_rad, __node.canopy_rt, __node.soil_opt, __node.wl_set, __node.rt_con, __node.rt_dim);
    _node_1x.can_rad = deepcopy(__node.can_rad);

    return [SIF_740(_node.can_rad,_node.wl_set) for _node in nodes]
end
