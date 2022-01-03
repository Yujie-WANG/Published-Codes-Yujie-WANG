###############################################################################
#
# Update g_sw using different stomatal models
#
# This part has been moved to StomataModels. Yet, keep this because the
# StomataModels version also addresses nighttime transpiration
#
###############################################################################
"""
Update stomatal conductance (`gsw`) for a given time period. This function
    simulates non-steady state `gsw`. However, `update_gsw!` is looped for 10+
    times to give steady state `gsw` value. `update_gsw!` is abstractized to
    use with both empirical and optimization stomatal models.
"""
function update_gsw! end




"""
This method is meant to work for any empirical stomatal model. Yet, in current
    project, we only tested the Ball-Berry model and Medlyn model. Use of this
    method using other stomatal models need to be cautious.

Note that we use a tune factor to downregulated photosynthetic capacitance as
    the β function to produce stomatal responses to drought, and thus the β
    factor is applied out of this function. As a result, a `1` was used at the
    placeholder of β in function `stomatal_conductance`.

    update_gsw!(clayer::CanopyLayer{FT},
                sm::EmpiricalStomatalModel{FT},
                photo_set::AbstractPhotoModelParaSet{FT},
                envir::AirLayer{FT},
                Δt::FT;
                τ::FT = FT(600),
                smoothing::Bool = false
    ) where {FT<:AbstractFloat}

Update stomatal conductance `g_sw` prognostically, given
- `clayer` A `CanopyLayer` type struct
- `sm` `EmpiricalStomatalModel` type stomatal model
- `photo_set` AbstractPhotoModelParaSet type photosynthesis model
- `envir` `AirLayer` type environmental conditions
- `Δt` Time interval for prognostic stomatal conductance
- `β` Beta factor for empirical models
- `τ` Time constant
- `smoothing` If true, use colimitation to smooth the An

We used a time constant of prognostic stomatal conductance of 600 s.
"""
update_gsw!(clayer::CanopyLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            Δt::FT;
            β::FT = FT(1),
            τ::FT = FT(600),
            smoothing::Bool = false
) where {FT<:AbstractFloat} =
(
    # unpack values
    @unpack g_bc, g_bw, g_lc, g_lw, g_m, g_sc, g_sw, n_leaf = clayer;

    # calculate steady state values
    # assume τ = 10 minutes
    for _iLF in 1:n_leaf
        # TODO add the bug fix to SPAC or StomataModels
        # bug fix: upadte APAR per leaf
        clayer.ps.APAR = clayer.APAR[_iLF];
        leaf_ETR!(photo_set, clayer.ps);

        _gsw_ss = max(sm.g0, stomatal_conductance(sm, clayer, envir, β, _iLF));
        g_sw[_iLF] += (_gsw_ss - g_sw[_iLF]) / τ * Δt;

        # bug fix: update g_lw and others as well
        g_lw[_iLF] = 1 / ( 1/g_sw[_iLF]  + 1/g_bw[_iLF] );
        g_sc[_iLF] = g_sw[_iLF] / FT(1.6);
        g_lc[_iLF] = 1 / ( 1/g_sc[_iLF] + 1/g_m[_iLF] + 1/g_bc[_iLF] );
    end;

    return nothing
)




"""
This method works for Wang et al. (2020) model only:


    update_gsw!(clayer::CanopyLayer{FT},
                sm::OSMWang{FT},
                photo_set::AbstractPhotoModelParaSet{FT},
                envir::AirLayer{FT},
                Δt::FT;
                β::FT = FT(600),
                τ::FT = FT(1e-6),
                smoothing::Bool = false
    ) where {FT<:AbstractFloat}

Update stomatal conductance `g_sw` prognostically, given
- `clayer` A `CanopyLayer` type struct
- `sm` `OSMWang` type stomatal model
- `photo_set` AbstractPhotoModelParaSet type photosynthesis model
- `envir` `AirLayer` type environmental conditions
- `Δt` Time interval for prognostic stomatal conductance
- `β` Beta factor for empirical models
- `τ` Time constant
- `smoothing` If true, use colimitation to smooth the An

We assumed that `Δgsw = (∂A∂E - ∂Θ∂E) * Δt * 1e-6` for OSMWang nodel.
"""
update_gsw!(clayer::CanopyLayer{FT},
            sm::OSMWang{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            Δt::FT;
            β::FT = FT(1),
            τ::FT = FT(1e-6),
            smoothing::Bool = false
) where {FT<:AbstractFloat} =
(
    # unpack values
    @unpack An, ec, g_bc, g_bw, g_lc, g_lw, g_m, g_max, g_min, g_sc, g_sw, n_leaf, p_sat, ps = clayer;
    @unpack p_atm, p_H₂O = envir;

    # update g_sw
    for iLF in 1:n_leaf
        # TODO add the bug fix to SPAC or StomataModels
        # bug fix: upadte APAR per leaf
        clayer.ps.APAR = clayer.APAR[iLF];
        leaf_ETR!(photo_set, clayer.ps);

        _gsw = g_sw[iLF] .+ FT(0.001);
        _glw = 1 / ( 1/_gsw + 1/g_bw[iLF] );
        _gsc = _gsw / FT(1.6);
        _glc = 1 / ( 1/_gsc + 1/g_m[iLF] + 1/g_bc[iLF] );
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), _glc);
        # TODO add colimitation back to Photosynthesis
        if smoothing
            _a1   = ps.Ac;
            _a2   = ps.Aj;
            _a3   = ps.Ap;
            _ag   = lower_quadratic(0.99, -(_a1 + _a2), _a1 * _a2);
            _ag   = lower_quadratic(0.99, -(_ag + _a3), _ag * _a3);
            ps.Ag = _ag;
            ps.An = ps.Ag - ps.Rd;
        end;
        _∂A   = ps.An - An[iLF];
        _e0   = g_lw[iLF] * (p_sat - p_H₂O) / p_atm;
        _e1   = _glw * (p_sat - p_H₂O) / p_atm;
        _∂E   = _e1 - _e0;
        _∂A∂E = _∂A / _∂E;
        _∂Θ∂E = FT(max(0.1, An[iLF]) / max(ec - _e0, eps(FT)));

        # disable nocturnal transpiration
        if clayer.ps.APAR < 1
            _∂Θ∂E = FT(Inf);
        end;

        # ensure that dgsw does not change too rapidly
        _Δgsw = (_∂A∂E - _∂Θ∂E) * τ * Δt;
        if _Δgsw > 0
            _Δgsw = min(_Δgsw, (g_max-g_sw[iLF]) / 4);
        else
            _Δgsw = max(_Δgsw, (g_min-g_sw[iLF]) / 4);
        end;
        g_sw[iLF] += _Δgsw;

        # bug fix: update g_lw and others as well
        g_lw[iLF] = 1 / ( 1/g_sw[iLF]  + 1/g_bw[iLF] );
        g_sc[iLF] = g_sw[iLF] / FT(1.6);
        g_lc[iLF] = 1 / ( 1/g_sc[iLF] + 1/g_m[iLF] + 1/g_bc[iLF] );
    end;

    return nothing
)








###############################################################################
#
# Simulate the diurnal cycle for a given set of environmental conditions
#
# Parts of this function has been included in the SoilPlantAirContinuum.jl
#     release. Yet, keep this function for now because nighttime transpiration
#     is not yet considerred in this project
#
###############################################################################
"""
Simulate the carbon, water, and energy fluxes.
"""
function simulation! end




"""
This method simulates the carbon, water, and energy fluxes for different modes
    include IJKX (original RT mode), 2KX (sunlit and shaded fraction per layer
    mode), KX (per layer mode), 2X (sunlit and shaded fraction mode), and 1X
    (global average mode).

    simulation!(proj::CanopyComplexity2021{FT},
                nodes::Vector{SPACMono{FT}},
                envir::AirLayer{FT};
                τ::FT = FT(1e-6),
                debugging::Bool = false
    ) where {FT<:AbstractFloat}

Simulate the gas exchange for different RT modes, given
- `proj` [`CanopyComplexity2021`](@ref) type project identifier
- `nodes` Different nodes
- `envir` Environmental conditions
- `τ` Time constant for the prognostic stomatal model
- `debugging` If true, display the debugging information
"""
function simulation!(
            proj::CanopyComplexity2021{FT},
            nodes::Vector{SPACMono{FT}},
            envir::AirLayer{FT};
            τ::FT = FT(1e-6),
            debugging::Bool = false
) where {FT<:AbstractFloat}
    # update environmental coonditions
    _node_ijkx, _node_2kx, _node_kx, _node_2x, _node_1x = nodes;
    for _node in nodes
        for _i in eachindex(_node.envirs)
            _node.envirs[_i] = deepcopy(envir);
        end
    end

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
    for _i_can in 1:_node_ijkx.n_canopy
        _iRT    = _node_ijkx.n_canopy + 1 - _i_can;
        _f_view = (can_opt.Ps[_iRT] + can_opt.Ps[_iRT+1]) / 2;

        # IJKX mode
        _iPS = _node_ijkx.plant_ps[_i_can];
        for _iLF in 1:_nSL
            _iPS.APAR[_iLF] = can_rad.absPAR_sunCab[(_iRT-1)*_nSL+_iLF] * FT(1e6);
            _iPS.LAIx[_iLF] = _f_view * f_SL[_iLF];
        end;
        _iPS.APAR[end] = can_rad.absPAR_shadeCab[_iRT] * FT(1e6);
        _iPS.LAIx[end] = 1 - _f_view;

        # 2KX mode
        _iPS         = _node_2kx.plant_ps[_i_can];
        _iPS.APAR[1] = mean(_node_ijkx.plant_ps[_i_can].APAR[1:end-1]);
        _iPS.LAIx[1] = _f_view;
        _iPS.APAR[2] = _node_ijkx.plant_ps[_i_can].APAR[end];
        _iPS.LAIx[2] = 1 - _f_view;

        # KX mode
        _iPS         = _node_kx.plant_ps[_i_can];
        _iPS.APAR[1] = _node_2kx.plant_ps[_i_can].APAR[1] * _f_view + _node_2kx.plant_ps[_i_can].APAR[2] * (1 - _f_view);
    end

    # update the PAR per layer for 2X and 1X mode
    _apar_sls = [_node_2kx.plant_ps[_i_can].APAR[1] for _i_can in 1:_node_2kx.n_canopy];
    _apar_shs = [_node_2kx.plant_ps[_i_can].APAR[2] for _i_can in 1:_node_2kx.n_canopy];
    _laix_sls = [_node_2kx.plant_ps[_i_can].LAIx[1] for _i_can in 1:_node_2kx.n_canopy];
    _laix_shs = [_node_2kx.plant_ps[_i_can].LAIx[2] for _i_can in 1:_node_2kx.n_canopy];
    _apar_all = [_node_kx.plant_ps[_i_can].APAR[1]  for _i_can in 1:_node_kx.n_canopy ];
    _node_2x.plant_ps[1].APAR[1]     = _apar_sls' * _laix_sls / sum(_laix_sls);
    _node_2x.plant_ps[2].APAR[1]     = _apar_shs' * _laix_shs / sum(_laix_shs);
    _node_2x.plant_ps[1].LAIx[1]     = 1;
    _node_2x.plant_ps[2].LAIx[1]     = 1;
    _node_1x.plant_ps[1].APAR[1]     = mean(_apar_all);
    _node_1x.plant_ps[1].LAIx[1]     = 1;
    _node_2x.plant_ps[1].LA          = _node_1x.plant_ps[1].LA * mean(_laix_sls);
    _node_2x.plant_ps[2].LA          = _node_1x.plant_ps[1].LA * mean(_laix_shs);
    _node_2x.plant_hs.leaves[1].area = _node_2x.plant_ps[1].LA;
    _node_2x.plant_hs.leaves[2].area = _node_2x.plant_ps[2].LA;

    # update the average Vcmax for 2X and 1X mode
    _vcmaxs = [_node_kx.plant_ps[_i_can].ps.Vcmax25 for _i_can in 1:_node_kx.n_canopy];
    _x_sl   = _laix_sls' * _vcmaxs / sum(_laix_sls) / _node_kx.plant_ps[end].ps.Vcmax25;
    _x_sh   = _laix_shs' * _vcmaxs / sum(_laix_shs) / _node_kx.plant_ps[end].ps.Vcmax25;
    _x_al   = mean(_vcmaxs) / _node_kx.plant_ps[end].ps.Vcmax25;
    _reps   = _node_ijkx.plant_ps[end].ps;
    _node_2x.plant_ps[1].ps.Vcmax25 = _x_sl * _reps.Vcmax25;
    _node_2x.plant_ps[1].ps.Vcmax   = _x_sl * _reps.Vcmax;
    _node_2x.plant_ps[1].ps.Jmax25  = _x_sl * _reps.Jmax25;
    _node_2x.plant_ps[1].ps.Jmax    = _x_sl * _reps.Jmax;
    _node_2x.plant_ps[1].ps.Rd25    = _x_sl * _reps.Rd25;
    _node_2x.plant_ps[1].ps.Rd      = _x_sl * _reps.Rd;
    _node_2x.plant_ps[2].ps.Vcmax25 = _x_sh * _reps.Vcmax25;
    _node_2x.plant_ps[2].ps.Vcmax   = _x_sh * _reps.Vcmax;
    _node_2x.plant_ps[2].ps.Jmax25  = _x_sh * _reps.Jmax25;
    _node_2x.plant_ps[2].ps.Jmax    = _x_sh * _reps.Jmax;
    _node_2x.plant_ps[2].ps.Rd25    = _x_sh * _reps.Rd25;
    _node_2x.plant_ps[2].ps.Rd      = _x_sh * _reps.Rd;
    _node_1x.plant_ps[1].ps.Vcmax25 = _x_al * _reps.Vcmax25;
    _node_1x.plant_ps[1].ps.Vcmax   = _x_al * _reps.Vcmax;
    _node_1x.plant_ps[1].ps.Jmax25  = _x_al * _reps.Jmax25;
    _node_1x.plant_ps[1].ps.Jmax    = _x_al * _reps.Jmax;
    _node_1x.plant_ps[1].ps.Rd25    = _x_al * _reps.Rd25;
    _node_1x.plant_ps[1].ps.Rd      = _x_al * _reps.Rd;
    if debugging
        @show _node_1x.plant_ps[1].ps.Vcmax25;

        # verify that sum of APAR equals everywhere
        for _i_node in 1:3
            _node = nodes[_i_node];
            @show sum([_node.plant_ps[_i].APAR' * _node.plant_ps[_i].LAIx * _node.plant_ps[_i].LAI for _i in 1:_node.n_canopy]);
            if _i_node > 1
                for _ps in _node.plant_ps
                    @show _ps.APAR;
                    @show _ps.LAIx;
                end
            end
        end
        for _i_node in 4:5
            _node = nodes[_i_node];
            @show sum(_node.plant_ps[1].APAR' * _node.plant_ps[1].LAIx * _node.plant_ps[1].LAI);
            for _ps in _node.plant_ps
                @show _ps.APAR;
                @show _ps.LAIx;
            end
        end
    end;

    # run gas exchange simulations for all nodes
    _ntimes = 720;
    _δt::FT = 14400 / _ntimes
    for _node in nodes
        for _ntime in 1:_ntimes
            _f_H₂O = 0;
            _f_CO₂ = 0;
            for _i_can in 1:_node.n_canopy
                _iEN = _node.envirs[_i_can];
                _iHS = _node.plant_hs.leaves[_i_can];
                _iPS = _node.plant_ps[_i_can];

                # update parameters such as critical flow
                update_leaf_TP!(_node.photo_set, _iPS, _iHS, _iEN);
                temperature_effects!(_iHS, (_iPS).T);

                # calculate the photosynthetic rates
                gas_exchange!(_node.photo_set, _iPS, _iEN, GswDrive());
                for _iLF in eachindex((_iPS).An)
                    _a1 = (_iPS).Ac[_iLF];
                    _a2 = (_iPS).Aj[_iLF];
                    _a3 = (_iPS).Ap[_iLF];
                    _ag = lower_quadratic(0.99, -(_a1 + _a2), _a1 * _a2);
                    _ag = lower_quadratic(0.99, -(_ag + _a3), _ag * _a3);
                    (_iPS).Ag[_iLF] = _ag;
                    (_iPS).An[_iLF] = _ag - (_iPS).ps.Rd;
                end;
                update_gsw!(_iPS, _node.stomata_model, _node.photo_set, _iEN, _δt; τ=τ, smoothing=true);
                gsw_control!(_node.photo_set, _iPS, _iEN);

                # update the flow rates
                _f_CO₂ += (_iPS).An' * (_iPS).LAIx * (_iPS).LA;
                _f_H₂O += (_iPS).g_lw' * (_iPS).LAIx * (_iPS).LA * ((_iPS).p_sat - (_iEN).p_H₂O) / _iEN.p_atm;
            end
            _node.f_npp = _f_CO₂ / _node.ga;
            _node.f_H₂O = _f_H₂O / _node.ga;

            # update flow profile and pressure history along the tree
            for _i_can in 1:_node.n_canopy
                _iEN      = _node.envirs[_i_can];
                _iPS      = _node.plant_ps[_i_can];
                _iLF      = _node.plant_hs.leaves[_i_can];
                _iLF.flow = _iPS.g_lw' * _iPS.LAIx * (_iPS.p_sat - _iEN.p_H₂O) / _iEN.p_atm;
            end
            for _iRT in _node.plant_hs.roots
                _iRT.flow = _f_H₂O / length(_node.plant_hs.roots);
            end

            # do not update history for now
            pressure_profile!(_node.plant_hs, SteadyStateMode(); update=false);

            # update canopy layer p_ups, which will be passed to each leaf
            for _i_can in 1:_node.n_canopy
                _iHS       = _node.plant_hs.leaves[_i_can];
                _iPS       = _node.plant_ps[_i_can];
                _iPS.p_ups = _iHS.p_ups;
            end
            # @show _ntime, [_iPS.p_ups for _iPS in _node.plant_ps];
        end
    end

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

    if debugging
        println("\n");
        @show _node_1x.plant_ps[1].ps.φs;
        @show _node_1x.plant_ps[1].ps.p_i;
        @show _node_1x.plant_ps[1].ps.Ac;
        @show _node_1x.plant_ps[1].ps.Aj;
        @show _node_1x.plant_ps[1].ps.Ap;
        @show _node_1x.plant_ps[1].ps.An;
        @show _node_1x.plant_ps[1].ps.Ag;
        @show _node_1x.plant_ps[1].ps.J;
        @show _node_1x.plant_ps[1].ps.J_pot;
        @show _node_1x.plant_ps[1].ps.Ja;
        @show _node_1x.plant_ps[1].ps.APAR;
        println("\n");
    end;

    return nothing
end
