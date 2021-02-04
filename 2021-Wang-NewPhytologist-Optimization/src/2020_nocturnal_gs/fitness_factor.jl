###############################################################################
#
# Calculate dAdE numerically by (1) find the Vcmax and (2) calculate dAdE
#
###############################################################################
"""
    calculate_dAdE!(
                leaf::Leaf{FT},
                env::AirLayer{FT},
                psm::AbstractPhotoModelParaSet{FT},
                a_n::FT,
                c_a::FT,
                e_leaf::FT,
                gbw::FT,
                p_atm::FT,
                par::FT,
                t_leaf::FT,
                vp_air::FT,
                dTdE::Bool
    ) where {FT<:AbstractFloat}

Calculate the marginal water use efficiency by firstly fitting the
    photosynthetic capacity and secondly numerically computing dAdE, given
- `leaf` Leaf type struct that stores leaf photosynthetic parameters
- `env` AirLayer type struct that stores environmental conditions
- `psm` AbstractPhotoModelParaSet type struct that stores photosynthesis
    model parameters, like the temperature dependencies
- `a_n` Given net photosynthetic rate
- `c_a` Given CO₂ concentration in ppm
- `e_leaf` Given leaf transpiration rate
- `gbw` Given leaf boundary layer conductance for H₂O
- `p_atm` Given atmospheric pressure
- `par` Given photosynthetically active radiation
- `t_leaf` Given leaf temperature
- `vp_air` Given vapor pressure in the air
- `dTdE::Bool` If true, account for leaf cooling from transpiration
"""
function calculate_dAdE!(
            leaf::Leaf{FT},
            env::AirLayer{FT},
            psm::AbstractPhotoModelParaSet{FT},
            a_n::FT,
            c_a::FT,
            e_leaf::FT,
            gbw::FT,
            p_atm::FT,
            par::FT,
            t_leaf::FT,
            vp_air::FT,
            dTdE::Bool
) where {FT<:AbstractFloat}
    env.p_atm = p_atm * 1000;
    env.p_H₂O = vp_air;
    env.p_a   = c_a * FT(1e-6) * p_atm * 1000;
    leaf.APAR   = par;
    leaf.T      = t_leaf + FT(273.15);
    leaf_temperature_dependence!(psm, leaf, env);

    tlf = leaf.T;
    glw = e_leaf / (leaf.p_sat - env.p_H₂O) * env.p_atm;
    gsw = 1 / (1/glw - 1/gbw);
    gsc = gsw / FT(1.6);
    gbc = gbw / FT(1.35);
    glc = 1 / (1/gsc + 1/gbc);

    # fit photosynthetic capacity first
    count = 0;
    while true
        count += 1;
        leaf_temperature_dependence!(psm, leaf, env);
        leaf_photo_from_glc!(psm, leaf, env, glc);

        if (abs(leaf.An - a_n) < 1e-5) || count>50
            break
        else
            leaf.Vcmax25 /= leaf.An / a_n;
            leaf.Jmax25  /= leaf.An / a_n;
            leaf.Rd25    /= leaf.An / a_n;
        end
    end

    # temperature related
    de    ::FT = FT(1e-5);
    lambda::FT = latent_heat_vapor(tlf) / 1000 * 18;
    cp    ::FT = FT(29.3);
    gbe   ::FT = 0.189 * sqrt(3.0/(0.72*0.1));
    fview ::FT = 0.2;
    ϵ     ::FT = 0.97;
    σ     ::FT = 5.67e-8;
    dTdE_low = 2 * cp * gbe + 4*fview*σ*ϵ*tlf^3;

    # comment this to disable temperature change due to de
    if dTdE
        tlf -= lambda / dTdE_low * de;
    end

    leaf_temperature_dependence!(psm, leaf, env, tlf);

    # calculate the dAdE
    a   = leaf.An;
    e   = e_leaf;
    ede = e + de;
    glw = ede / (leaf.p_sat - env.p_H₂O) * env.p_atm;
    gsw = 1 / (1/glw - 1/gbw);
    gsc = gsw / FT(1.6);
    gbc = gbw / FT(1.35);
    glc = 1 / (1/gsc + 1/gbc);
    leaf_photo_from_glc!(psm, leaf, env, glc);
    a_de = leaf.An;

    dAdE = (a_de - a) / de;

    return dAdE
end








###############################################################################
#
# Calculate dRdE analytically
#
###############################################################################
"""
    calculate_dRdE(
                psm::AbstractPhotoModelParaSet{FT},
                r_n::FT,
                t_leaf::FT,
                wind::FT
    ) where {FT<:AbstractFloat}

Calculate marginal respiratory reduction, given
- `psm` AbstractPhotoModelParaSet type struct that stores photosynthesis
    model parameters, like the temperature dependencies
- `r_n` Given leaf respiration rate
- `t_leaf` Given leaf temperature
- `wind` Given wind speed
"""
function calculate_dRdE(
            psm::AbstractPhotoModelParaSet{FT},
            r_n::FT,
            t_leaf::FT,
            wind::FT
) where {FT<:AbstractFloat}
    tlf  = t_leaf + FT(273.15);
    dRdE = calculate_dTdE(t_leaf, wind) * r_n * psm.ΓsT.ΔHa_to_R / tlf^2;

    return dRdE
end




function calculate_dTdE(
            t_leaf::FT,
            wind::FT
) where {FT<:AbstractFloat}
    tlf = t_leaf + FT(273.15);

    lambda = latent_heat_vapor(tlf) / 1000 * 18;
    cp     = FT(29.3);
    gbe    = FT(0.189 * sqrt(wind/(0.72*0.1)));
    fview  = FT(1 / 4.7585);
    emis   = FT(0.97);
    σ      = FT(5.67e-8);
    dRdE_low = 2 * cp * gbe + 4 * fview * σ * emis * tlf^3;

    dTdE = lambda / dRdE_low;

    return dTdE
end




function calculate_Tn(
            psm::AbstractPhotoModelParaSet{FT},
            r_25::FT,
            drde::FT,
            wind::FT
) where {FT<:AbstractFloat}
    # x in Celsius
    @inline f(x) = (
        r_n = r_25 * temperature_correction(psm.ReT, x+FT(273.15));
        tar = calculate_dRdE(psm, r_n, x, wind);
        return -abs(tar - drde);
    );

    ms  = BisectionMethod{FT}(x_min=15, x_max=35);
    st  = SolutionTolerance{FT}(1e-3, 50);
    sol = find_peak(f, ms, st);

    return sol
end




function calculate_Gn(
            psm::AbstractPhotoModelParaSet{FT},
            r_25::FT,
            drde::FT,
            wind::FT,
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    t_leaf = calculate_Tn(psm, r_25, drde, wind);
    cp     = FT(29.3);
    gbe    = FT(0.189 * sqrt(wind/(0.72*0.1)));
    lambda = latent_heat_vapor(t_leaf + FT(273.15)) / 1000 * 18;
    e_leaf = -2 * cp * gbe * (t_leaf+FT(273.15) - envir.t_air) / lambda;
    g_leaf = e_leaf / (saturation_vapor_pressure(t_leaf+FT(273.15)) - envir.p_H₂O) * envir.p_atm;

    return g_leaf
end








###############################################################################
#
# Calculate fitness factor numerically
#
###############################################################################
"""
    fitness_factor(
                project::NocturnalGS2020{FT};
                dTdE::Bool,
                display_result::Bool
    ) where {FT<:AbstractFloat}

Calculate the mean and standard deviation of the fitness factors, given
- `project` NocturnalGS2020 type project control
- `dTdE` Whether accounting for leaf cooling effect for daytime photosynthesis
- `display_result` Optional. If true, display the results
"""
function fitness_factor(
            project::NocturnalGS2020{FT};
            dTdE::Bool = true,
            display_result::Bool = true
) where {FT<:AbstractFloat}
    # create leaf and environment
    leaf = Leaf{FT}();
    env  = AirLayer{FT}();
    psm  = C3CLM(FT);

    # type in the data
    list_an = FT[14.13694195, 13.42925102, 12.87522973,
                 11.22623193, 9.703102473, 8.929087186];
    list_ca = FT[378.9937429, 379.8834581, 380.5374590,
                 382.6031277, 384.1396535, 385.3239880];
    list_el = FT[0.009035248, 0.008752303, 0.008786093,
                 0.008478174, 0.009116404, 0.008558623];
    list_gb = FT[2.042710619, 2.042936103, 2.043917173,
                 2.043629211, 2.042936793, 2.043330535];
    list_pa = FT[86.35421905, 86.35338095, 86.35479524,
                 86.35650952, 86.35798095, 86.35594286];
    list_pr = FT[600.1000000, 600.0120000, 599.9910000,
                 600.0590000, 600.0560000, 599.7760000];
    list_tl = FT[31.32444762, 31.82125714, 32.91324286,
                 33.78967619, 33.14503333, 33.84225238];
    list_vp = FT[2.547731134, 2.692166152, 2.790227162,
                 2.807653129, 2.880491612, 2.830240939];
    list_rn = FT[1.140798360, 1.003127783, 1.013440461,
                 1.058663616, 0.882447345, 0.824649337];
    list_tn = FT[27.50960443, 27.66397390, 27.79030838,
                 28.24151591, 28.31068597, 28.51245940];

    list_dade = calculate_dAdE!.([leaf],
                                 [env],
                                 [psm],
                                 list_an,
                                 list_ca,
                                 list_el,
                                 list_gb,
                                 list_pa,
                                 list_pr,
                                 list_tl,
                                 list_vp,
                                 [dTdE]);
    list_drde = calculate_dRdE.( [psm],
                                 list_rn,
                                 list_tn,
                                 FT[0.1]);
    list_ff = list_drde ./ list_dade;

    if display_result
        @show list_dade;
        @show list_drde;
        @show list_ff;
        @show mean(list_ff);
        @show std(list_ff);
    end

    return nothing
end
