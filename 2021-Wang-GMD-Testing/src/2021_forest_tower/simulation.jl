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
                τ::FT = FT(600)
    ) where {FT<:AbstractFloat}

Update stomatal conductance `g_sw` prognostically, given
- `clayer` A `CanopyLayer` type struct
- `sm` `EmpiricalStomatalModel` type stomatal model
- `photo_set` AbstractPhotoModelParaSet type photosynthesis model
- `envir` `AirLayer` type environmental conditions
- `Δt` Time interval for prognostic stomatal conductance
- `τ` Time constant

We used a time constant of prognostic stomatal conductance of 600 s.
"""
update_gsw!(clayer::CanopyLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            Δt::FT;
            τ::FT = FT(600)
) where {FT<:AbstractFloat} =
(
    # unpack values
    @unpack g_sw, n_leaf, ps = clayer;

    # calculate steady state values
    # assume τ = 10 minutes
    for _iLF in 1:n_leaf
        _gsw_ss = max(sm.g0,
                      stomatal_conductance(sm, clayer, envir, FT(1), _iLF));
        g_sw[_iLF] += (_gsw_ss - g_sw[_iLF]) / τ * Δt;
    end;

    return nothing
)




"""
This method works for Wang et al. (2020) model only:


    update_gsw!(clayer::CanopyLayer{FT},
                sm::OSMWang{FT},
                photo_set::AbstractPhotoModelParaSet{FT},
                envir::AirLayer{FT},
                Δt::FT
    ) where {FT<:AbstractFloat}

Update stomatal conductance `g_sw` prognostically, given
- `clayer` A `CanopyLayer` type struct
- `sm` `OSMWang` type stomatal model
- `photo_set` AbstractPhotoModelParaSet type photosynthesis model
- `envir` `AirLayer` type environmental conditions
- `Δt` Time interval for prognostic stomatal conductance

We assumed that `Δgsw = (∂A∂E - ∂Θ∂E) * Δt * 1e-6` for OSMWang nodel.
"""
update_gsw!(clayer::CanopyLayer{FT},
            sm::OSMWang{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            Δt::FT;
            τ::FT = FT(1e-6)
) where {FT<:AbstractFloat} =
(
    # unpack values
    @unpack An, ec, g_bc, g_bw, g_lw, g_m, g_max, g_min, g_sw, n_leaf, p_sat,
            ps = clayer;
    @unpack p_atm, p_H₂O = envir;

    # update g_sw
    for iLF in 1:n_leaf
        _gsw = g_sw[iLF] .+ FT(0.001);
        _glw = 1 / ( 1/_gsw + 1/g_bw[iLF] );
        _gsc = _gsw / FT(1.6);
        _glc = 1 / ( 1/_gsc + 1/g_m[iLF] + 1/g_bc[iLF] );
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), _glc);
        _∂A   = ps.An - An[iLF];
        _e0   = g_lw[iLF] * (p_sat - p_H₂O) / p_atm;
        _e1   = _glw * (p_sat - p_H₂O) / p_atm;
        _∂E   = _e1 - _e0;
        _∂A∂E = _∂A / _∂E;
        _∂Θ∂E = FT(max(0.1, An[iLF]) / max(ec - _e0, 1e-7));

        # ensure that dgsw does not change too rapidly
        _Δgsw = (_∂A∂E - _∂Θ∂E) * τ * Δt;
        if _Δgsw > 0
            _Δgsw = min(_Δgsw, (g_max-g_sw[iLF]) / 4);
        else
            _Δgsw = max(_Δgsw, (g_min-g_sw[iLF]) / 4);
        end;
        g_sw[iLF] += _Δgsw;
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
Simulate the carbon and water flux, and save the fluxes into the input data
    frame. The simulated data may be used to compare with the flux tower
    estimates, and thus be used to fit site level variables like photosynthetic
    capacity and hydraulic conductance.
"""
function simulation! end




"""
This method simulates the carbon and water fluxes for flux tower sites Niwot
    Ridge and Ozark. To reduce the number of uncertainty, some parameters are
    prescribed, such as soil moisture and leaf temperature.

    simulation!(proj::ForestTower2021{FT},
                node::SPACMono{FT},
                df::DataFrame,
                site::String,
                rbase::Q10TD{FT} = Q10TD{FT}(0, 298.15, 1.7)
    ) where {FT<:AbstractFloat}

Run diurnal cycle of land model simulation, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `node` `SPACMono` type struct from `SoilPlantAirContinuum` module
- `df` Dataframe of weather data. Note that simulated data are stored in `df`
- `site` Which flux tower site
- `rbase` Q10 based temperature dependency struct. Default is 0.

Note that zenith angle is computed at each time step, and the light environment
    is treated as 0 if the zenith angle is higher than 88°. Also, because there
    is no direct and diffuse light data available, the partition of the two is
    done by adjsting incomming PPFD to observed PPFD.
"""
simulation!(proj::ForestTower2021{FT},
            node::SPACMono{FT},
            df::DataFrame,
            site::String,
            rbase::Q10TD{FT} = Q10TD{FT}(0, 298.15, 1.7)
) where {FT<:AbstractFloat} =
(
    # 0.1 unpack data
    @unpack angles, can_opt, can_rad, canopy_rt, envirs, f_SL, ga, in_rad,
            latitude, leaves_rt, n_canopy, photo_set, plant_hs, plant_ps,
            rt_con, soil_opt, stomata_model, wl_set = node;
    @unpack lidf, nAzi, nIncl = canopy_rt;
    @unpack dWL, iPAR = wl_set;
    in_rad_bak = deepcopy(in_rad);
    nSL = nAzi * nIncl;
    in_Erad = in_rad_bak.E_direct .+ in_rad_bak.E_diffuse;
    in_PPFD = sum( e2phot(dWL, in_Erad)[iPAR] ) * FT(1e6);

    # iterate through the weather data
    for i in eachindex(df.Day)
        # update soil water matrices
        w_soil = FT(df.SWC[i]) / 100;
        t_soil = FT(df.T_SOIL[i] + 273.15);
        # need to adjust SWC to avoid problem in residual SWC at Niwot Ridge
        if site == "NiwotRidge"
            ψ_soil = soil_p_25_swc(plant_hs.roots[1].sh,
                                   w_soil + plant_hs.roots[1].sh.Θr);
        else
            ψ_soil = soil_p_25_swc(plant_hs.roots[1].sh, w_soil);
        end;
        for root in plant_hs.roots
            root.p_ups = ψ_soil;
        end;

        # update PAR related information
        zenith = zenith_angle(latitude, FT(df.Day[i]), FT(df.Hour[i]),
                              FT(df.Minu[i]));
        zenith = min(88, zenith);
        angles.sza = zenith;
        in_rad.E_direct  .= in_rad_bak.E_direct  .* df.PPFD[i] ./ in_PPFD;
        in_rad.E_diffuse .= in_rad_bak.E_diffuse .* df.PPFD[i] ./ in_PPFD;
        canopy_geometry!(canopy_rt, angles, can_opt, rt_con);
        canopy_matrices!(leaves_rt, can_opt);
        short_wave!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, rt_con);
        canopy_fluxes!(canopy_rt, can_opt, can_rad, in_rad, soil_opt,
                       leaves_rt, wl_set, rt_con);

        # update fluxes
        f_H₂O = 0;
        f_CO₂ = 0;
        for i_can in 1:n_canopy
            iEN = envirs[i_can];
            iHS = plant_hs.leaves[i_can];
            iPS = plant_ps[i_can];
            iRT = n_canopy + 1 - i_can;

            # update environmental conditions
            iEN.t_air = df.T_AIR[i] + 273.15;
            iEN.p_atm = df.P_ATM[i] * 1000;
            iEN.p_a   = iEN.p_atm * 4e-4;
            iEN.p_O₂  = iEN.p_atm * 0.209;
            iEN.p_sat = saturation_vapor_pressure(iEN.t_air);
            iEN.vpd   = df.VPD[i] * 100;
            iEN.p_H₂O = iEN.p_sat - iEN.vpd;
            iEN.RH    = iEN.p_H₂O / iEN.p_sat;
            iEN.wind  = df.WIND[i];

            # prescribe leaf temperature
            _tl = (df.LW_OUT[i] / 0.97 / K_STEFAN() ) ^ 0.25;
            iPS.T = _tl;
            update_leaf_TP!(photo_set, iPS, iHS, iEN);
            temperature_effects!(iHS, FT(_tl));

            # calculate the fraction of sunlit and shaded leaves
            f_view = (can_opt.Ps[iRT] + can_opt.Ps[iRT+1]) / 2;
            for iLF in 1:nSL
                iPS.APAR[iLF] = can_rad.absPAR_sunCab[(iRT-1)*nSL+iLF] *
                                FT(1e6);
                iPS.LAIx[iLF] = f_view * f_SL[iLF];
            end;
            iPS.APAR[end] = can_rad.absPAR_shadeCab[iRT] * FT(1e6);
            iPS.LAIx[end] = 1 - f_view;

            # iterate for 15 times to find steady state solution
            for iter in 1:15
                # calculate the photosynthetic rates
                gas_exchange!(photo_set, iPS, iEN, GswDrive());
                update_gsw!(iPS, stomata_model, photo_set, iEN, FT(120));
                gsw_control!(photo_set, iPS, iEN);
            end;

            # update the flow rates
            for iLF in 1:(nSL+1)
                f_CO₂ += iPS.An[iLF] * iPS.LAIx[iLF] * iPS.LA;
                f_H₂O += iPS.g_lw[iLF] * (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm *
                         iPS.LAIx[iLF] * iPS.LA;
            end;
        end;

        # update flow profile and pressure history along the tree
        for i_can in 1:n_canopy
            iEN = envirs[i_can];
            iLF = plant_hs.leaves[i_can];
            iPS = plant_ps[i_can];
            iST = plant_hs.branch[i_can];
            iLF.flow = sum(iPS.g_lw .* iPS.LAIx) *
                       (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm;
            iST.flow = iLF.flow * iPS.LA;
        end;
        plant_hs.trunk.flow = sum([iST.flow for iST in plant_hs.branch]);
        for iRT in plant_hs.roots
            iRT.flow = plant_hs.trunk.flow / length(plant_hs.roots);
        end;

        # do not update history for now
        pressure_profile!(plant_hs, SteadyStateMode(); update=false);

        # update canopy layer p_ups, which will be passed to each leaf
        for _i_can in 1:n_canopy
            _iHS = plant_hs.leaves[_i_can];
            _iPS = plant_ps[_i_can];
            _iPS.p_ups = _iHS.p_ups;
        end

        # update Vcmax based on leaf water potential for empirical models
        if !(typeof(stomata_model) <: OSMWang{FT})
            for i_can in 1:n_canopy
                iEN = envirs[i_can];
                iHS = plant_hs.leaves[i_can];
                iPS = plant_ps[i_can];
                iRT = n_canopy + 1 - i_can;

                _ratio = xylem_risk(iHS, iHS.flow);
                update_VJR!(node, _ratio);
                iPS.T_old = 0;
            end;
        end;

        # calculate respiration from other components
        _r = photo_TD_from_set(rbase, t_soil);

        # update the data frame
        df.ModC[i] = f_CO₂ / ga - _r;
        df.ModE[i] = f_H₂O / ga;
    end;

    return nothing
)




"""
This method update the stomatal conductance to steady state so as to plot the
    stomatal responses to different environmental stimuli:

    simulation!(proj::ForestTower2021{FT},
                node::SPACMono{FT},
                beta::Bool
    ) where {FT<:AbstractFloat}

Simulate steady state stomatal conductance, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `node` `SPACMono` type struct from `SoilPlantAirContinuum` module
- `beta` If true, photosynthetic capacity is updated for empirical model
"""
simulation!(proj::ForestTower2021{FT},
            node::SPACMono{FT},
            beta::Bool
) where {FT<:AbstractFloat} =
(
    # 0.1 unpack data
    @unpack angles, can_opt, can_rad, canopy_rt, envirs, f_SL, ga, in_rad,
            latitude, leaves_rt, n_canopy, photo_set, plant_hs, plant_ps,
            rt_con, soil_opt, stomata_model, wl_set = node;
    @unpack lidf, nAzi, nIncl = canopy_rt;
    @unpack dWL, iPAR = wl_set;
    nSL = nAzi * nIncl;

    # update PAR related information
    canopy_geometry!(canopy_rt, angles, can_opt, rt_con);
    canopy_matrices!(leaves_rt, can_opt);
    short_wave!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, rt_con);
    canopy_fluxes!(canopy_rt, can_opt, can_rad, in_rad, soil_opt,
                   leaves_rt, wl_set, rt_con);

    # update fluxes
    f_H₂O = 0;
    f_CO₂ = 0;
    for i_can in 1:n_canopy
        iEN = envirs[i_can];
        iHS = plant_hs.leaves[i_can];
        iPS = plant_ps[i_can];
        iRT = n_canopy + 1 - i_can;

        # prescribe leaf temperature
        update_leaf_TP!(photo_set, iPS, iHS, iEN);
        temperature_effects!(iHS, iPS.T);

        # calculate the fraction of sunlit and shaded leaves
        f_view = (can_opt.Ps[iRT] + can_opt.Ps[iRT+1]) / 2;
        for iLF in 1:nSL
            iPS.APAR[iLF] = can_rad.absPAR_sunCab[(iRT-1)*nSL+iLF] *
                            FT(1e6);
            iPS.LAIx[iLF] = f_view * f_SL[iLF];
        end;
        iPS.APAR[end] = can_rad.absPAR_shadeCab[iRT] * FT(1e6);
        iPS.LAIx[end] = 1 - f_view;

        # iterate for 60 times (60 s interval) to find steady state solution
        for iter in 1:100
            # calculate the photosynthetic rates
            gas_exchange!(photo_set, iPS, iEN, GswDrive());
            if !(typeof(stomata_model) <: OSMWang{FT})
                update_gsw!(iPS, stomata_model, photo_set, iEN, FT(60));
            else
                update_gsw!(iPS, stomata_model, photo_set, iEN, FT(60);
                            τ=FT(1e-7));
            end
            gsw_control!(photo_set, iPS, iEN);

            # update flow profile and pressure history along the tree
            for i_can in 1:n_canopy
                iEN = envirs[i_can];
                iLF = plant_hs.leaves[i_can];
                iPS = plant_ps[i_can];
                iST = plant_hs.branch[i_can];
                iLF.flow = sum(iPS.g_lw .* iPS.LAIx) *
                           (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm;
                iST.flow = iLF.flow * iPS.LA;
            end;
            plant_hs.trunk.flow = sum([iST.flow for iST in plant_hs.branch]);
            for iRT in plant_hs.roots
                iRT.flow = plant_hs.trunk.flow / length(plant_hs.roots);
            end;

            # do not update history for now
            pressure_profile!(plant_hs, SteadyStateMode(); update=false);

            # update canopy layer p_ups, which will be passed to each leaf
            for _i_can in 1:n_canopy
                _iHS = plant_hs.leaves[_i_can];
                _iPS = plant_ps[_i_can];
                _iPS.p_ups = _iHS.p_ups;
            end;

            # update Vcmax based on leaf water potential for empirical models
            if !(typeof(stomata_model) <: OSMWang{FT}) && beta
                for i_can in 1:n_canopy
                    iEN = envirs[i_can];
                    iHS = plant_hs.leaves[i_can];
                    iPS = plant_ps[i_can];
                    iRT = n_canopy + 1 - i_can;

                    _ratio = xylem_risk(iHS, iHS.flow);
                    update_VJR!(node, _ratio);
                    iPS.T_old = 0;
                end;
            end;
        end;

        # update the flow rates
        for iLF in 1:(nSL+1)
            f_CO₂ += iPS.An[iLF] * iPS.LAIx[iLF] * iPS.LA;
            f_H₂O += iPS.g_lw[iLF] * (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm *
                     iPS.LAIx[iLF] * iPS.LA;
        end;
    end;

    node.f_npp = f_CO₂;
    node.f_H₂O = f_H₂O;

    return nothing
)




"""
Compared to the method above, this method has an extra parameter in the
    function input variables: `ft_patm`.

    simulation!(proj::ForestTower2021{FT},
                node::SPACMono{FT},
                df::DataFrame,
                site::String,
                ft_patm::Bool,
                rbase::Q10TD{FT} = Q10TD{FT}(0, 298.15, 1.7)
    ) where {FT<:AbstractFloat}

Run diurnal cycle of land model simulation, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `node` `SPACMono` type struct from `SoilPlantAirContinuum` module
- `df` Dataframe of weather data. Note that simulated data are stored in `df`
- `site` Which flux tower site
- `ft_patm` If true, use Flux tower `p_atm`; else, use sea level `p_atm`
- `rbase` Q10 based temperature dependency struct. Default is 0.

Given that this function differs only in `p_atm`, we run the simulation by
    changing the `p_atm` in the data frame to sea level mean.
"""
simulation!(proj::ForestTower2021{FT},
            node::SPACMono{FT},
            df::DataFrame,
            site::String,
            ft_patm::Bool,
            rbase::Q10TD{FT} = Q10TD{FT}(0, 298.15, 1.7)
) where {FT<:AbstractFloat} =
(
    # if true, change the p_atm in the data frame to 101.325 kPa
    if ft_patm
        df.P_ATM .= 101.325;
    end;
    simulation!(proj, node, df, site, rbase);

    return nothing
)








###############################################################################
#
# Calculate monthly mean SIF
#
###############################################################################
"""
Simulate the time series of SIF. Note that as this simulation is meant to
    compare with TROPOMI SIF retrievals, we simulate the SIF at satallite
    overpassing time only, namely from 12:45 PM to 14:30 PM.
"""
function sif_series! end




"""
This method is modified from the [`simulation!`](!ref). The differences are
    that we simulated the SIF only at the satallite overpassing time and that
    we only save SIF in the data frame but neglect carbon and water fluxes.

    sif_series!(proj::ForestTower2021{FT},
                node::SPACMono{FT},
                df::DataFrame,
                year::Int,
                site::String,
                sNPQ::Bool
    ) where {FT<:AbstractFloat}

Run annual SIF simulations, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `node` `SPACMono` type struct from `SoilPlantAirContinuum` module
- `df` Dataframe of weather data. Note that simulated data are stored in `df`
- `year` Which year to simulate
- `sNPQ` If true, use prescribed sustained NPQ; otherwise, use 0. Note that
    this `true` option only works for Niwot Ridge
"""
sif_series!(proj::ForestTower2021{FT},
            node::SPACMono{FT},
            df::DataFrame,
            year::Int,
            site::String,
            sNPQ::Bool
) where {FT<:AbstractFloat} =
(
    # recorded median sustained NPQ, data from Raczka et al. (2019)
    if sNPQ
        @assert site == "NiwotRidge";
        npqs = FT[5.28, 6.32, 6.36, 4.74, 1.70, 0.12, 0.02, 0.02, 0.04, 0.54,
                  2.44, 5.36];
    else
        npqs = zeros(FT,12);
    end;

    # 0.1 unpack data
    @unpack angles, can_opt, can_rad, canopy_rt, envirs,f_SL,  ga, in_rad,
            latitude, leaves_rt, n_canopy, photo_set, plant_hs, plant_ps,
            rt_con, rt_dim, soil_opt, stomata_model, wl_set = node;
    @unpack lidf, nAzi, nIncl = canopy_rt;
    @unpack dWL, iPAR = wl_set;
    in_rad_bak = deepcopy(in_rad);
    nSL = nAzi * nIncl;
    in_Erad = in_rad_bak.E_direct .+ in_rad_bak.E_diffuse;
    in_PPFD = sum( e2phot(dWL, in_Erad)[iPAR] ) * FT(1e6);
    angles.vza = 0;

    df.Month  = [parse(Int, string(ts)[5:6]) for ts in df.TIME];
    df.SIF740 = [NaN for j in eachindex(df.TIME)];
    df.susNPQ = [NaN for j in eachindex(df.TIME)];
    for mon in 1:12
        mask = (df.Month .== mon);
        df.susNPQ[mask] .= npqs[mon];
    end;

    # iterate through the weather data
    @info tinfo("Simulating SIF for $(site) in $(year)...");
    @showprogress for i in eachindex(df.Day)
        # TROPOMI SIF overpass time is from 12:45PM to 2:30 PM
        # Thus, we choose data within this time window only
        _hour = df.Hour[i] + df.Minu[i]/60;
        _calc = (13 < _hour < 14);

        # run the simulation only if _calc is true
        if _calc
            # update soil water matrices
            w_soil = FT(df.SWC[i]) / 100;
            t_soil = FT(df.T_SOIL[i] + 273.15);
            # need to adjust SWC to avoid problem in residual SWC
            if site == "NiwotRidge"
                ψ_soil = soil_p_25_swc(plant_hs.roots[1].sh,
                                       w_soil + plant_hs.roots[1].sh.Θr);
            else
                ψ_soil = soil_p_25_swc(plant_hs.roots[1].sh, w_soil);
            end;
            for root in plant_hs.roots
                root.p_ups = ψ_soil;
            end;

            # update relative azimuth angle from time
            angles.raa = 180 - (df.Hour[i] + df.Minu[i]/60) * 15;

            # update PAR related information
            zenith = zenith_angle(latitude, FT(df.Day[i]), FT(df.Hour[i]),
                                  FT(df.Minu[i]));
            zenith = min(88, zenith);
            angles.sza = zenith;
            in_rad.E_direct  .= in_rad_bak.E_direct  .* df.PPFD[i] ./ in_PPFD;
            in_rad.E_diffuse .= in_rad_bak.E_diffuse .* df.PPFD[i] ./ in_PPFD;
            canopy_geometry!(canopy_rt, angles, can_opt, rt_con);
            canopy_matrices!(leaves_rt, can_opt);
            short_wave!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, rt_con);
            canopy_fluxes!(canopy_rt, can_opt, can_rad, in_rad, soil_opt,
                           leaves_rt, wl_set, rt_con);

            # update fluxes
            for i_can in 1:n_canopy
                iEN = envirs[i_can];
                iHS = plant_hs.leaves[i_can];
                iPS = plant_ps[i_can];
                iRT = n_canopy + 1 - i_can;

                # update sustained NPQ
                iPS.ps.Ks = df.susNPQ[i];

                # update environmental conditions
                iEN.t_air = df.T_AIR[i] + 273.15;
                iEN.p_atm = df.P_ATM[i] * 1000;
                iEN.p_a   = iEN.p_atm * 4e-4;
                iEN.p_O₂  = iEN.p_atm * 0.209;
                iEN.p_sat = saturation_vapor_pressure(iEN.t_air);
                iEN.vpd   = df.VPD[i] * 100;
                iEN.p_H₂O = iEN.p_sat - iEN.vpd;
                iEN.RH    = iEN.p_H₂O / iEN.p_sat;
                iEN.wind  = df.WIND[i];

                # prescribe leaf temperature
                _tl = (df.LW_OUT[i] / 0.97 / K_STEFAN()) ^ 0.25;
                iPS.T = _tl;
                update_leaf_TP!(photo_set, iPS, iHS, iEN);
                temperature_effects!(iHS, FT(_tl));

                # calculate the fraction of sunlit and shaded leaves
                f_view = (can_opt.Ps[iRT] + can_opt.Ps[iRT+1]) / 2;
                for iLF in 1:nSL
                    iPS.APAR[iLF] = can_rad.absPAR_sunCab[(iRT-1)*nSL+iLF] *
                                    FT(1e6);
                    iPS.LAIx[iLF] = f_view * f_SL[iLF];
                end;
                iPS.APAR[end] = can_rad.absPAR_shadeCab[iRT] * FT(1e6);
                iPS.LAIx[end] = 1 - f_view;

                # iterate for 15 times to find steady state solution
                for iter in 1:15
                    # calculate the photosynthetic rates
                    gas_exchange!(photo_set, iPS, iEN, GswDrive());
                    update_gsw!(iPS, stomata_model, photo_set, iEN, FT(120));
                    gsw_control!(photo_set, iPS, iEN);
                end;

                # update fluorescence quantum yield from leaf photosynthesis
                can_rad.ϕ_sun[:,:,iRT] .= reshape(view(iPS.ϕs,1:nSL), nIncl,
                                                  nAzi);
                can_rad.ϕ_shade[iRT] = iPS.ϕs[end];
            end;

            # do SIF simulation
            SIF_fluxes!(leaves_rt, can_opt, can_rad, canopy_rt, soil_opt,
                        wl_set, rt_con, rt_dim);

            # update SIF from the table if zenith angle <= 70
            if zenith <= 70
                df.SIF740[i] = SIF_740(can_rad, wl_set);
            end;

            # update flow profile and pressure history along the tree
            for i_can in 1:n_canopy
                iEN = envirs[i_can];
                iLF = plant_hs.leaves[i_can];
                iPS = plant_ps[i_can];
                iST = plant_hs.branch[i_can];
                iLF.flow = sum(iPS.g_lw .* iPS.LAIx) *
                           (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm;
                iST.flow = iLF.flow * iPS.LA;
            end;
            plant_hs.trunk.flow = sum([iST.flow for iST in plant_hs.branch]);
            for iRT in plant_hs.roots
                iRT.flow = plant_hs.trunk.flow / length(plant_hs.roots);
            end;
            pressure_profile!(plant_hs, SteadyStateMode(); update=true);

            # update canopy layer p_ups, which will be passed to each leaf
            for _i_can in 1:n_canopy
                _iHS = plant_hs.leaves[_i_can];
                _iPS = plant_ps[_i_can];
                _iPS.p_ups = _iHS.p_ups;
            end

            # update Vcmax based on leaf water potential for empirical models
            if !(typeof(stomata_model) <: OSMWang{FT})
                for i_can in 1:n_canopy
                    iEN = envirs[i_can];
                    iHS = plant_hs.leaves[i_can];
                    iPS = plant_ps[i_can];
                    iRT = n_canopy + 1 - i_can;

                    _ratio = xylem_risk(iHS, iHS.flow);
                    update_VJR!(node, _ratio);
                    iPS.T_old = 0;
                end;
            end;
        end;
    end;

    # [nanmean( df.SIF740[df.Month .== mon]) for mon in 1:12]
    _folder = "/home/wyujie/RAID/Data/FLUXNET2015/simulation";
    if sNPQ
        save_csv!("$(_folder)/sif_$(site)_$(year)_sNPQ.csv", df);
    else
        save_csv!("$(_folder)/sif_$(site)_$(year).csv", df);
    end;

    return nothing
)




"""
This method is the case when we simulate SIF to match per TROPOMI observation.
    Note that the data is reprocessed data using [`match_TROPOMI`](@ref):

    sif_series!(proj::ForestTower2021{FT},
                node::SPACMono{FT},
                df::DataFrame,
                year::Int,
                site::String;
                modis_lai::Bool = false
    ) where {FT<:AbstractFloat}

Simulate SIF time series to match each TROPOMI SIF observation, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `node` `SPACMono` type struct from `SoilPlantAirContinuum` module
- `df` Dataframe of weather data. Note that simulated data are stored in `df`
- `year` Which year to simulate
- `modis_lai` Optional. If true, use MODIS LAI time series
"""
sif_series!(proj::ForestTower2021{FT},
            node::SPACMono{FT},
            df::DataFrame,
            year::Int,
            site::String;
            modis_lai::Bool = false
) where {FT<:AbstractFloat} =
(
    # 0.1 unpack data
    @unpack angles, can_opt, can_rad, canopy_rt, envirs,f_SL,  ga, in_rad,
            latitude, leaves_rt, n_canopy, photo_set, plant_hs, plant_ps,
            rt_con, rt_dim, soil_opt, stomata_model, wl_set = node;
    @unpack lidf, nAzi, nIncl = canopy_rt;
    @unpack dWL, iPAR = wl_set;
    in_rad_bak = deepcopy(in_rad);
    nSL = nAzi * nIncl;
    in_Erad = in_rad_bak.E_direct .+ in_rad_bak.E_diffuse;
    in_PPFD = sum( e2phot(dWL, in_Erad)[iPAR] ) * FT(1e6);
    angles.vza = 0;

    df.Month  = [parse(Int, string(ts)[5:6]) for ts in df.TIME];
    df.SIF740 = [NaN for j in eachindex(df.TIME)];

    # iterate through the weather data
    @info tinfo("Simulating SIF for $(site) in $(year)...");
    @showprogress for i in eachindex(df.Day)
        # update LAI
        if modis_lai update_LAI!(node, FT(df.LAI[i])); end;

        # update soil water matrices
        w_soil = FT(df.SWC[i]) / 100;
        t_soil = FT(df.T_SOIL[i] + 273.15);
        # need to adjust SWC to avoid problem in residual SWC
        if site == "NiwotRidge"
            ψ_soil = soil_p_25_swc(plant_hs.roots[1].sh,
                                   w_soil + plant_hs.roots[1].sh.Θr);
        else
            ψ_soil = soil_p_25_swc(plant_hs.roots[1].sh, w_soil);
        end;
        for root in plant_hs.roots
            root.p_ups = ψ_soil;
        end;

        # update relative azimuth angle from time

        # update PAR related information
        angles.raa = df.RAA[i];
        angles.sza = df.SZA[i];
        angles.vza = df.VZA[i];
        in_rad.E_direct  .= in_rad_bak.E_direct  .* df.PPFD[i] ./ in_PPFD;
        in_rad.E_diffuse .= in_rad_bak.E_diffuse .* df.PPFD[i] ./ in_PPFD;
        canopy_geometry!(canopy_rt, angles, can_opt, rt_con);
        canopy_matrices!(leaves_rt, can_opt);
        short_wave!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, rt_con);
        canopy_fluxes!(canopy_rt, can_opt, can_rad, in_rad, soil_opt,
                       leaves_rt, wl_set, rt_con);

        # update fluxes
        for i_can in 1:n_canopy
            iEN = envirs[i_can];
            iHS = plant_hs.leaves[i_can];
            iPS = plant_ps[i_can];
            iRT = n_canopy + 1 - i_can;

            # update environmental conditions
            iEN.t_air = df.T_AIR[i] + 273.15;
            iEN.p_atm = df.P_ATM[i] * 1000;
            iEN.p_a   = iEN.p_atm * 4e-4;
            iEN.p_O₂  = iEN.p_atm * 0.209;
            iEN.p_sat = saturation_vapor_pressure(iEN.t_air);
            iEN.vpd   = df.VPD[i] * 100;
            iEN.p_H₂O = iEN.p_sat - iEN.vpd;
            iEN.RH    = iEN.p_H₂O / iEN.p_sat;
            iEN.wind  = df.WIND[i];

            # prescribe leaf temperature
            _tl = (df.LW_OUT[i] / 0.97 / K_STEFAN()) ^ 0.25;
            iPS.T = _tl;
            update_leaf_TP!(photo_set, iPS, iHS, iEN);
            temperature_effects!(iHS, FT(_tl));

            # calculate the fraction of sunlit and shaded leaves
            f_view = (can_opt.Ps[iRT] + can_opt.Ps[iRT+1]) / 2;
            for iLF in 1:nSL
                iPS.APAR[iLF] = can_rad.absPAR_sunCab[(iRT-1)*nSL+iLF] *
                                FT(1e6);
                iPS.LAIx[iLF] = f_view * f_SL[iLF];
            end;
            iPS.APAR[end] = can_rad.absPAR_shadeCab[iRT] * FT(1e6);
            iPS.LAIx[end] = 1 - f_view;

            # iterate for 15 times to find steady state solution
            for iter in 1:15
                # calculate the photosynthetic rates
                gas_exchange!(photo_set, iPS, iEN, GswDrive());
                update_gsw!(iPS, stomata_model, photo_set, iEN, FT(120));
                gsw_control!(photo_set, iPS, iEN);
            end;

            # update fluorescence quantum yield from leaf photosynthesis
            can_rad.ϕ_sun[:,:,iRT] .= reshape(view(iPS.ϕs,1:nSL), nIncl,
                                              nAzi);
            can_rad.ϕ_shade[iRT] = iPS.ϕs[end];
        end;

        # do SIF simulation
        SIF_fluxes!(leaves_rt, can_opt, can_rad, canopy_rt, soil_opt, wl_set,
                    rt_con, rt_dim);

        # update SIF from the table
        df.SIF740[i] = SIF_740(can_rad, wl_set);

        # update flow profile and pressure history along the tree
        for i_can in 1:n_canopy
            iEN = envirs[i_can];
            iLF = plant_hs.leaves[i_can];
            iPS = plant_ps[i_can];
            iST = plant_hs.branch[i_can];
            iLF.flow = sum(iPS.g_lw .* iPS.LAIx) *
                       (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm;
            iST.flow = iLF.flow * iPS.LA;
        end;
        plant_hs.trunk.flow = sum([iST.flow for iST in plant_hs.branch]);
        for iRT in plant_hs.roots
            iRT.flow = plant_hs.trunk.flow / length(plant_hs.roots);
        end;
        pressure_profile!(plant_hs, SteadyStateMode(); update=true);

        # update canopy layer p_ups, which will be passed to each leaf
        for _i_can in 1:n_canopy
            _iHS = plant_hs.leaves[_i_can];
            _iPS = plant_ps[_i_can];
            _iPS.p_ups = _iHS.p_ups;
        end

        # update Vcmax based on leaf water potential for empirical models
        if !(typeof(stomata_model) <: OSMWang{FT})
            for i_can in 1:n_canopy
                iEN = envirs[i_can];
                iHS = plant_hs.leaves[i_can];
                iPS = plant_ps[i_can];
                iRT = n_canopy + 1 - i_can;

                _ratio = xylem_risk(iHS, iHS.flow);
                update_VJR!(node, _ratio);
                iPS.T_old = 0;
            end;
        end;
    end;

    # [nanmean( df.SIF740[df.Month .== mon]) for mon in 1:12]
    _folder = "/home/wyujie/RAID/Data/FLUXNET2015/simulation";
    if modis_lai
        save_csv!("$(_folder)/sif_$(site)_$(year)_LAI.csv", df);
    else
        save_csv!("$(_folder)/sif_$(site)_$(year).csv", df);
    end;

    return nothing
)




"""
This function is a wrapper to run four SIF simulations, and they are:
- Site: Niwot Ridge, Year: 2018
- Site: Niwot Ridge, Year: 2019
- Site: Ozark, Year: 2018
- Site: Ozark, Year: 2019

    sif_series!(proj::ForestTower2021{FT}) where {FT<:AbstractFloat}

Run annual SIF simulations, given
- `proj` [`ForestTower2021`](@ref) type project identifier
"""
sif_series!(proj::ForestTower2021{FT}) where {FT<:AbstractFloat} =
(
    # read data from the fitting results
    _folder = "/home/wyujie/RAID/Data/FLUXNET2015/simulation";
    _fitting = read_csv("$(_folder)/fitted_all.csv");

    # simulate model SIF for Niwot Ridge
    _df_18  = match_TROPOMI(proj, "NiwotRidge", 2018);
    _df_19  = match_TROPOMI(proj, "NiwotRidge", 2019);
    _m_18   = (_fitting.Site .== "NiwotRidge") .* (_fitting.Model .== "osm") .*
              (_fitting.Year .== 2018);
    _m_19   = (_fitting.Site .== "NiwotRidge") .* (_fitting.Model .== "osm") .*
              (_fitting.Year .== 2019);
    _osm_18 = create_spac(proj, "NiwotRidge", OSMWang{FT}(),
                          FT(_fitting.Vcmax[_m_18][1]),
                          FT(_fitting.K[_m_18][1]) / 2,
                          FT(46.82));
    _osm_19 = create_spac(proj, "NiwotRidge", OSMWang{FT}(),
                          FT(_fitting.Vcmax[_m_19][1]),
                          FT(_fitting.K[_m_19][1]) / 2,
                          FT(46.82));

    # run sif simulations
    sif_series!(proj, deepcopy(_osm_18), _df_18, 2018, "NiwotRidge");
    sif_series!(proj, deepcopy(_osm_19), _df_19, 2019, "NiwotRidge");
    sif_series!(proj, deepcopy(_osm_18), _df_18, 2018, "NiwotRidge";
                modis_lai=true);
    sif_series!(proj, deepcopy(_osm_19), _df_19, 2019, "NiwotRidge";
                modis_lai=true);

    # simulate model SIF for Ozark
    _df_18  = match_TROPOMI(proj, "Ozark", 2018);
    _df_18  = match_TROPOMI(proj, "Ozark", 2019);
    _m_18   = (_fitting.Site .== "Ozark") .* (_fitting.Model .== "osm") .*
              (_fitting.Year .== 2018);
    _m_19   = (_fitting.Site .== "Ozark") .* (_fitting.Model .== "osm") .*
              (_fitting.Year .== 2019);
    _osm_18 = create_spac(proj, "Ozark", OSMWang{FT}(),
                          FT(_fitting.Vcmax[_m_18][1]),
                          FT(_fitting.K[_m_18][1]) / 2,
                          FT(57.23));
    _osm_19 = create_spac(proj, "Ozark", OSMWang{FT}(),
                          FT(_fitting.Vcmax[_m_19][1]),
                          FT(_fitting.K[_m_19][1]) / 2,
                          FT(57.23));

    # run sif simulations
    sif_series!(proj, deepcopy(_osm_18), _df_18, 2018, "Ozark");
    sif_series!(proj, deepcopy(_osm_19), _df_19, 2019, "Ozark");
    sif_series!(proj, deepcopy(_osm_18), _df_18, 2018, "Ozark";
                modis_lai=true);
    sif_series!(proj, deepcopy(_osm_19), _df_19, 2019, "Ozark";
                modis_lai=true);

    return nothing;
)
