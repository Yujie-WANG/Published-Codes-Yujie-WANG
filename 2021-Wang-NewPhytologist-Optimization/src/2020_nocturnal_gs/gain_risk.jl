###############################################################################
#
# Optimize nighttime flow
#
###############################################################################
function optimal_en_gn(
            proj::NocturnalGS2020{FT},
            spac::SPACSimple{FT},
            can::CanopyLayer{FT},
            psm::AbstractPhotoModelParaSet{FT},
            dsm::OSMWang{FT},
            ff::FT = FT(0.1);
            ΔT::FT = FT(0)
) where {FT<:AbstractFloat}
    # function to find zero
    @inline f(x) = (
        gr = marginal_gain_risk!(proj, spac, can, psm, dsm, x; ΔT=ΔT);
        return gr[1] - gr[2]*ff;
    );

    # solve for optimum
    ec = critical_flow(spac.hs) / spac.laba;
    ms = NewtonBisectionMethod{FT}(1e-6, ec * 0.9, 1e-6);
    st = SolutionTolerance{FT}(1e-6, 50);
    en = find_zero(f, ms, st);
    gn = en / (can.ps.p_sat - spac.envir.p_H₂O) * spac.envir.p_atm;

    return en, gn
end








###############################################################################
#
# Calculate marginal gain and risk
#
###############################################################################
function marginal_gain_risk!(
            proj::NocturnalGS2020{FT},
            spac::SPACSimple{FT},
            can::CanopyLayer{FT},
            psm::AbstractPhotoModelParaSet{FT},
            dsm::OSMWang{FT},
            flow::FT;
            ΔT::FT = FT(0)
) where {FT<:AbstractFloat}
    # calculate temperature at nighttime
    can.T = leaf_temperature(spac, FT(0), flow*spac.laba);

    # use nighttime temperature to calculate A and R
    leaf_photo_from_envir!(psm, can, spac.hs, spac.envir, dsm);

    # calculate marginal gain and risk
    ∂R∂E = dRdE(proj, spac, can, psm);

    # update temperatures
    if ΔT != 0
        can.T            += ΔT;
        spac.envir.t_air += ΔT;
        spac.envir.p_sat  = saturation_vapor_pressure(spac.envir.t_air);
        spac.envir.vpd    = spac.envir.p_sat - spac.envir.p_H₂O;
        leaf_photo_from_envir!(psm, can, spac.hs, spac.envir, dsm);
    end
    ∂Θ∂E = dΘdE(proj, spac, can, psm, flow);

    # change the values back
    if ΔT != 0
        can.T            -= ΔT;
        spac.envir.t_air -= ΔT;
        spac.envir.p_sat  = saturation_vapor_pressure(spac.envir.t_air);
        spac.envir.vpd    = spac.envir.p_sat - spac.envir.p_H₂O;
    end

    return ∂R∂E, ∂Θ∂E
end




function dRdE(
            proj::NocturnalGS2020{FT},
            spac::SPACSimple{FT},
            can::CanopyLayer{FT},
            psm::AbstractPhotoModelParaSet{FT}
) where {FT<:AbstractFloat}
    @unpack ps, T = can;

    return dTdE(proj, spac, T) * ps.Rd * psm.ΓsT.ΔHa_to_R / T^2
end




function dTdE(
            proj::NocturnalGS2020{FT},
            spac::SPACSimple{FT},
            t_leaf::FT
) where {FT<:AbstractFloat}
    lambda = latent_heat_vapor(t_leaf) * MOLMASS_WATER(FT);
    cp     = FT(29.3);
    gbe    = FT(0.189 * sqrt(spac.envir.wind/(0.72*spac.width)));
    emis   = FT(0.97);
    denom  = 2 * cp * gbe + 4 / spac.lai * K_STEFAN(FT) * emis * t_leaf^3;
    ∂T∂E   = lambda / denom;

    return ∂T∂E
end




function dΘdE(
            proj::NocturnalGS2020{FT},
            spac::SPACSimple{FT},
            can::CanopyLayer{FT},
            psm::AbstractPhotoModelParaSet{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack envir = spac;
    @unpack ps = can;

    # calculate g_lc and a_net
    g_lw = flow / max(1, ps.p_sat - envir.p_H₂O) * envir.p_atm;
    g_bw = ps.g_bc * FT(1.35);
    g_sw = min(spac.g_max, 1 / max(FT(1e-3), 1/g_lw - 1/g_bw));
    g_sc = g_sw / FT(1.6);
    g_lc = 1 / (1/g_sc + 1/ps.g_bc);
    leaf_photo_from_glc!(psm, ps, envir, g_lc);

    # calculate e_crit (out of the function to speed up later)
    can.ec = critical_flow(spac.hs, spac.ec) / spac.laba;

    return ps.An / (can.ec - flow)
end
