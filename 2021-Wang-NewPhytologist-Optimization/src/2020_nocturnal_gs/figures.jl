###############################################################################
#
# Figure 1: model framework
#
###############################################################################
"""
    plot_model_framework(
                proj::NocturnalGS2020{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot model framework, given
- `proj` NocturnalGS2020 type project control
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and san-serif font
"""
function plot_model_framework(
            proj::NocturnalGS2020{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    # use sans and latex
    if use_latex use_serif_tex(); end

    # define variables to work with
    spac = create_spac(proj);
    can = CanopyLayer{FT}(n_leaf=1);
    psm = C3CLM(FT);
    dsm = OSMWang{FT}();
    can.ps.APAR = 1000;
    spac.envir.wind = 0.1;
    psm.ReT = RespirationTDBernacchi(FT);

    # generate arrays of gain and risk
    list_e = collect(FT, 1:50:2001) * FT(1e-6);
    list_a = similar(list_e);
    list_c = similar(list_e);
    list_r = similar(list_e);
    list_t = similar(list_e);
    list_b = similar(list_e);
    for i in eachindex(list_e)
        list_b[i],list_c[i] = marginal_gain_risk!(proj, spac, can, psm, dsm,
                                                  list_e[i]);
        list_a[i] = can.ps.An;
        list_r[i] = can.ps.Rd;
        list_t[i] = can.ps.T;
    end

    # plot the figure
    fig,axs = create_canvas("GN-1", ncol=3);
    ax1,ax2,ax3 = axs;

    ax1.plot(list_e .* 1000, list_t .- 273.15, "k-");
    ax2.plot(list_e .* 1000, list_r, "r-", label=LS_Rn);
    ax2.plot(list_e .* 1000, list_a, "c-", label=LS_Ad);
    ax3.plot(list_e .* 1000, list_b, "r-", label=LS_n_∂Rn∂En);
    ax3.plot(list_e .* 1000, list_c ./ 10, "c-", label=LS_ff_∂Θd∂Ed);
    ax3.plot(1.0764, 144.94, "ko", markersize=10);
    ax2.legend(loc="center right");
    ax3.legend(loc="lower right");

    set_xlims!(axs, [[0,2] for i in 1:3]);
    set_ylims!(axs, [[22.9,24], [0,6.6], [0,240]]);
    set_xlabels!(axs, [LS_E_unit_m for i in 1:3]);
    set_ylabels!(axs, [LS_Tleaf_unit,
                       LS_Rn * " or " * LS_Ad * "\n" * latex_unit("A"),
                       "Margian gain or cost\n" * latex_unit("WUE")]);
    set_titles!(axs, loc="left");

    # save figure
    fig.set_tight_layout(true);
    if saving
        fig.savefig("figures/2020_nocturnal_gs/1_model_framework.pdf",
                bbox_inches="tight");
    end

    return nothing
end








###############################################################################
#
# Figure 2: model prediction
#
###############################################################################
"""
    plot_model_prediction(
                proj::NocturnalGS2020{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot model prediction, given
- `proj` NocturnalGS2020 type project control
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and san-serif font
"""
function plot_model_prediction(
            proj::NocturnalGS2020{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    # use sans and latex
    if use_latex use_serif_tex(); end

    # define variables to work with
    spac = create_spac(proj);
    can  = CanopyLayer{FT}(n_leaf=1);
    psm  = C3CLM(FT);
    dsm  = OSMWang{FT}();
    psm.ReT = RespirationTDBernacchi(FT);

    @inline initialize_vars() = (
        spac = create_spac(proj);
        can  = CanopyLayer{FT}(n_leaf=1);
        can.ps.APAR     = 1000;
        spac.envir.wind = 0.1;
    );

    # generate arrays of gain and risk
    list_e = collect(FT, 1:50:2201) * FT(1e-6);
    list_c = similar(list_e);
    list_b = similar(list_e);
    @inline generate_curves() = (
        for i in eachindex(list_e)
            list_b[i],list_c[i] = marginal_gain_risk!(proj, spac, can, psm,
                                                      dsm, list_e[i]);
        end
    );

    # plot the figure
    fig,axs = create_canvas("GN-2"; nrow=2, ncol=3);
    ax1,ax2,ax3,ax4,ax5,ax6 = axs;

    # plot baseline for all panels
    initialize_vars();
    generate_curves();
    for ax in axs
        ax.plot(list_e .* 1000, list_b, "r-");
        ax.plot(list_e .* 1000, list_c ./ 10, "c-");
    end

    # soil becomes drier
    initialize_vars();
    spac.p_soil        = -0.5;
    spac.hs.root.p_ups = spac.p_soil;
    generate_curves();
    ax1.plot(list_e .* 1000, list_c ./ 10, "c--", label="Drier soil");
    ax1.plot(1.076, 144.95, "o", color="gray" , markersize=8);
    ax1.plot(0.516, 145.20, "o", color="black", markersize=8);
    ax1.annotate("", xy=(0.54,152), xytext=(1.05,151), arrowprops=ARROW_K);

    # atmospheric CO₂ increase
    initialize_vars();
    spac.envir.p_a = 120;
    generate_curves();
    ax2.plot(list_e .* 1000, list_c ./ 10, "c--", label="Higher " * LS_Ca);
    ax2.plot(1.076, 144.95, "o", color="gray" , markersize=8);
    ax2.plot(0.247, 148.50, "o", color="black", markersize=8);
    ax2.annotate("", xy=(0.274,154), xytext=(1.05,151), arrowprops=ARROW_K);

    # atmospheric VPD increases
    initialize_vars();
    spac.envir.p_H₂O = 200;
    spac.envir.vpd   = spac.envir.p_sat - spac.envir.p_H₂O;
    generate_curves();
    ax3.plot(list_e .* 1000, list_c ./ 10, "c--", label="Higher VPD");
    ax3.plot(1.076, 144.95, "o", color="gray" , markersize=8);
    ax3.plot(1.354, 144.60, "o", color="black", markersize=8);
    ax3.annotate("", xy=(1.339,151), xytext=(1.09,151), arrowprops=ARROW_K);

    # respiration rate increases
    initialize_vars();
    can.ps.Rd25 = 1.3;
    generate_curves();
    ax4.plot(list_e .* 1000, list_b, "r--");
    ax4.plot(list_e .* 1000, list_c ./ 10, "c--", label="Higher " * LS_Rn);
    ax4.plot(1.076, 144.95, "o", color="gray" , markersize=8);
    ax4.plot(1.846, 184.54, "o", color="black", markersize=8);
    ax4.annotate("", xy=(1.747,188), xytext=(1.05,151), arrowprops=ARROW_K);

    # air temperature increases
    initialize_vars();
    spac.envir.t_air = 303.15;
    spac.envir.p_sat = saturation_vapor_pressure(spac.envir.t_air);
    spac.envir.p_H₂O = spac.envir.p_sat - spac.envir.vpd;
    generate_curves();
    ax5.plot(list_e .* 1000, list_b, "r--");
    ax5.plot(list_e .* 1000, list_c ./ 10, "c--", label="Higher T");
    ax5.plot(1.076, 144.95, "o", color="gray" , markersize=8);
    ax5.plot(2.122, 184.30, "o", color="black", markersize=8);
    ax5.annotate("", xy=(2.04,184), xytext=(1.18,151), arrowprops=ARROW_K);

    # fitness factor decreases
    initialize_vars();
    generate_curves();
    ax6.plot(list_e .* 1000, list_c .* 0.08, "c--", label="Lower " * LS_ff);
    ax6.plot(1.076, 144.95, "o", color="gray" , markersize=8);
    ax6.plot(1.640, 142.76, "o", color="black", markersize=8);
    ax6.annotate("", xy=(1.6,148), xytext=(1.15,151), arrowprops=ARROW_K);

    # add legend for each axis
    for ax in axs
        ax.legend(loc="lower right");
    end

    set_xlims!(axs, [0,2.2]);
    set_ylims!(axs, [0,240]);
    set_xlabels!([ax4,ax5,ax6], LS_E_unit_m);
    set_ylabels!([ax1,ax4], "Margian gain or cost\n" * latex_unit("WUE"));
    set_titles!(axs, loc="left");

    # save figure
    fig.set_tight_layout(true);
    if saving
        fig.savefig("figures/2020_nocturnal_gs/2_model_prediction.pdf",
                bbox_inches="tight");
    end

    return nothing
end








###############################################################################
#
# Figure 3-5: molde comparison with experimental observations
#
###############################################################################
"""
    plot_model_comparison(
                proj::NocturnalGS2020{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot model comparison with experimental observations, given
- `proj` NocturnalGS2020 type project control
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and san-serif font
"""
function plot_model_comparison(
            proj::NocturnalGS2020{FT};
            saving::Bool = false,
            use_latex::Bool=true
) where {FT<:AbstractFloat}
    # use serif and latex
    if use_latex use_serif_tex(); end

    #
    #
    # For Psoil, Ca, and VPD
    #
    #
    # read the data
    file_c = joinpath(artifact"2020_nocturnal_gs", "g_co2.csv");
    file_d = joinpath(artifact"2020_nocturnal_gs", "g_vpd.csv");
    file_p = joinpath(artifact"2020_nocturnal_gs", "g_ps.csv" );
    data_c = DataFrame(CSV.File(file_c));
    data_d = DataFrame(CSV.File(file_d));
    data_p = DataFrame(CSV.File(file_p));

    # define variables to work with
    spac = create_spac(proj);
    can  = CanopyLayer{FT}(n_leaf=1);
    psm  = C3CLM(FT);
    dsm  = OSMWang{FT}();
    psm.ReT = RespirationTDBernacchi(FT);

    @inline initialize_vars() = (
        spac = create_spac(proj);
        can  = CanopyLayer{FT}(n_leaf=1);
        can.ps.APAR     = 1000;
        spac.envir.wind = 0.1;
    );

    # plot the figure
    fig,axs = create_canvas("GN-3", ncol=3);
    ax1,ax2,ax3 = axs;

    # vs Psoil
    list_p = (data_p).Ps * FT(-0.1);
    list_g = (data_p).Gn;
    ax1.plot(list_p[ 1:  4], list_g[ 1:  4], color="k",
             mfc="none", marker="o", linestyle=":");
    ax1.plot(list_p[ 5:  8], list_g[ 5:  8], color="y",
             mfc="none", marker="o", linestyle=":");
    ax1.plot(list_p[ 9: 12], list_g[ 9: 12], color="g",
             mfc="none", marker="o", linestyle=":");
    ax1.plot(list_p[13: 16], list_g[13: 16], color="b",
             mfc="none", marker="o", linestyle=":");
    ax1.plot(list_p[17: 20], list_g[17: 20], color="c",
             mfc="none", marker="o", linestyle=":");
    ax1.plot(list_p[21:end], list_g[21:end], color="r",
             mfc="none", marker="o", linestyle=":");

    # genrate theoretical responses
    initialize_vars();
    mod_ps = collect(FT, 0:-0.1:-1.5);
    mod_gs = similar(mod_ps);
    for i in eachindex(mod_ps)
        spac.p_soil        = mod_ps[i];
        spac.hs.root.p_ups = mod_ps[i];
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, FT(0.15));
        mod_gs[i] = gn;
    end
    ax1.plot(mod_ps, mod_gs, color="silver", linewidth=10, zorder=1);

    # vs Ca
    list_c = (data_c).Ca;
    list_g = (data_c).Gn;
    ax2.plot(list_c[ 1:  5], list_g[ 1:  5], color="r",
             mfc="none", marker="o", linestyle=":");
    ax2.plot(list_c[ 6: 10], list_g[ 6: 10], color="y",
             mfc="none", marker="o", linestyle=":");
    ax2.plot(list_c[11: 15], list_g[11: 15], color="g",
             mfc="none", marker="o", linestyle=":");
    ax2.plot(list_c[16: 21], list_g[16: 21], color="b",
             mfc="none", marker="o", linestyle=":");
    ax2.plot(list_c[22: 27], list_g[22: 27], color="c",
             mfc="none", marker="o", linestyle=":");
    ax2.plot(list_c[28:end], list_g[28:end], color="k",
             mfc="none", marker="o", linestyle=":");

    # genrate theoretical responses
    initialize_vars();
    mod_cs = collect(FT, 5:5:80);
    mod_gs = similar(mod_cs);
    for i in eachindex(mod_cs)
        spac.envir.p_a = mod_cs[i];
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, FT(0.15));
        mod_gs[i] = gn;
    end
    ax2.plot(mod_cs .* 10, mod_gs, color="silver", linewidth=10, zorder=1);

    # vs D
    list_d = (data_d).D;
    list_g = (data_d).Gn;
    ax3.plot(list_d[ 1: 9], list_g[ 1: 9], color="r",
             mfc="none", marker="o", linestyle=":");
    ax3.plot(list_d[11:18], list_g[11:18], color="g",
             mfc="none", marker="o", linestyle=":");
    ax3.plot(list_d[20:28], list_g[20:28], color="b",
             mfc="none", marker="o", linestyle=":");
    ax3.plot(list_d[30:37], list_g[30:37], color="y",
             mfc="none", marker="o", linestyle=":");
    ax3.plot(list_d[39:47], list_g[39:47], color="c",
             mfc="none", marker="o", linestyle=":");
    ax3.plot(list_d[49:58], list_g[49:58], color="k",
             mfc="none", marker="o", linestyle=":");

    # genrate theoretical responses
    initialize_vars();
    mod_ds = collect(FT, 500:50:3000);
    mod_gs = similar(mod_ds);
    for i in eachindex(mod_ds)
        spac.envir.vpd   = mod_ds[i];
        spac.envir.p_H₂O = spac.envir.p_sat - spac.envir.vpd;
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, FT(0.15));
        mod_gs[i] = gn;
    end
    ax3.plot(mod_ds ./ 1000, mod_gs, color="silver", linewidth=10, zorder=1);

    # set axes ticks and labes
    set_ylims!(axs, [[0,0.08], [0,0.4], [0,0.13]]);
    set_xlabels!(axs, [LS_Psoil_unit, LS_Ca_unit, "VPD (kPa)"]);
    set_ylabels!([ax1], [LS_gwn_unit]);
    set_titles!(axs, loc="left");

    # save figure
    fig.set_tight_layout(true);
    if saving
        fig.savefig("figures/2020_nocturnal_gs/3_model_comparison.pdf",
                    bbox_inches="tight");
    end

    #
    #
    # Responses to Psoil, Ca, and VPD with ΔT
    #
    #
    # plot the figure
    fig,axs = create_canvas("GN-S2", ncol=3);
    ax1,ax2,ax3 = axs;

    # vs Psoil
    list_p = (data_p).Ps * FT(-0.1);
    list_g = (data_p).Gn;
    ax1.plot(list_p[ 1:  4], list_g[ 1:  4], color="k",
             mfc="none", marker="o", linestyle=":");
    ax1.plot(list_p[ 5:  8], list_g[ 5:  8], color="y",
             mfc="none", marker="o", linestyle=":");
    ax1.plot(list_p[ 9: 12], list_g[ 9: 12], color="g",
             mfc="none", marker="o", linestyle=":");
    ax1.plot(list_p[13: 16], list_g[13: 16], color="b",
             mfc="none", marker="o", linestyle=":");
    ax1.plot(list_p[17: 20], list_g[17: 20], color="c",
             mfc="none", marker="o", linestyle=":");
    ax1.plot(list_p[21:end], list_g[21:end], color="r",
             mfc="none", marker="o", linestyle=":");

    # genrate theoretical responses
    initialize_vars();
    mod_ps = collect(FT, 0:-0.1:-1.5);
    mod_gs = similar(mod_ps);
    for i in eachindex(mod_ps)
        spac.p_soil        = mod_ps[i];
        spac.hs.root.p_ups = mod_ps[i];
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, FT(0.15); ΔT=FT(10));
        mod_gs[i] = gn;
    end
    ax1.plot(mod_ps, mod_gs, color="silver", linewidth=10, zorder=1);

    # vs Ca
    list_c = (data_c).Ca;
    list_g = (data_c).Gn;
    ax2.plot(list_c[ 1:  5], list_g[ 1:  5], color="r",
             mfc="none", marker="o", linestyle=":");
    ax2.plot(list_c[ 6: 10], list_g[ 6: 10], color="y",
             mfc="none", marker="o", linestyle=":");
    ax2.plot(list_c[11: 15], list_g[11: 15], color="g",
             mfc="none", marker="o", linestyle=":");
    ax2.plot(list_c[16: 21], list_g[16: 21], color="b",
             mfc="none", marker="o", linestyle=":");
    ax2.plot(list_c[22: 27], list_g[22: 27], color="c",
             mfc="none", marker="o", linestyle=":");
    ax2.plot(list_c[28:end], list_g[28:end], color="k",
             mfc="none", marker="o", linestyle=":");

    # genrate theoretical responses
    initialize_vars();
    mod_cs = collect(FT, 10:5:80);
    mod_gs = similar(mod_cs);
    for i in eachindex(mod_cs)
        spac.envir.p_a = mod_cs[i];
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, FT(0.15); ΔT=FT(10));
        mod_gs[i] = gn;
    end
    ax2.plot(mod_cs .* 10, mod_gs, color="silver", linewidth=10, zorder=1);

    # vs D
    list_d = (data_d).D;
    list_g = (data_d).Gn;
    ax3.plot(list_d[ 1: 9], list_g[ 1: 9], color="r",
             mfc="none", marker="o", linestyle=":");
    ax3.plot(list_d[11:18], list_g[11:18], color="g",
             mfc="none", marker="o", linestyle=":");
    ax3.plot(list_d[20:28], list_g[20:28], color="b",
             mfc="none", marker="o", linestyle=":");
    ax3.plot(list_d[30:37], list_g[30:37], color="y",
             mfc="none", marker="o", linestyle=":");
    ax3.plot(list_d[39:47], list_g[39:47], color="c",
             mfc="none", marker="o", linestyle=":");
    ax3.plot(list_d[49:58], list_g[49:58], color="k",
             mfc="none", marker="o", linestyle=":");

    # genrate theoretical responses
    initialize_vars();
    mod_ds = collect(FT, 500:50:3000);
    mod_gs = similar(mod_ds);
    for i in eachindex(mod_ds)
        spac.envir.vpd   = mod_ds[i];
        spac.envir.p_H₂O = spac.envir.p_sat - spac.envir.vpd;
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, FT(0.15); ΔT=FT(10));
        mod_gs[i] = gn;
    end
    ax3.plot(mod_ds ./ 1000, mod_gs, color="silver", linewidth=10, zorder=1);

    # set axes ticks and labes
    set_ylims!(axs, [[0,0.08], [0,0.4], [0,0.13]]);
    set_xlabels!(axs, [LS_Psoil_unit, LS_Ca_unit, "VPD (kPa)"]);
    set_ylabels!([ax1], [LS_gwn_unit]);
    set_titles!(axs, loc="left");

    # save figure
    fig.set_tight_layout(true);
    if saving
        fig.savefig("figures/2020_nocturnal_gs/s2_model_comparison.pdf",
                    bbox_inches="tight");
    end

    #
    #
    # For temperature
    #
    #
    # read the data
    file_t = joinpath(artifact"2020_nocturnal_gs", "g_tem.csv");
    data_t = DataFrame(CSV.File(file_t));
    list_t = (data_t).Tl;
    list_g = (data_t).E ./ (data_t).D .* 86.0;
    list_f = list_g ./ (data_t).F;

    # plot the figure
    fig,axs = create_canvas("GN-4", ncol=2);
    ax1,ax2 = axs;

    # vs T not corrected
    ax1.plot(list_t[ 1: 18], list_g[ 1: 18], color="r",
             mfc="none", marker="o", linestyle="");
    ax1.plot(list_t[19: 35], list_g[19: 35], color="y",
             mfc="none", marker="o", linestyle="");
    ax1.plot(list_t[36: 58], list_g[36: 58], color="g",
             mfc="none", marker="o", linestyle="");
    ax1.plot(list_t[59: 74], list_g[59: 74], color="b",
             mfc="none", marker="o", linestyle="");
    ax1.plot(list_t[75: 90], list_g[75: 90], color="c",
             mfc="none", marker="o", linestyle="");
    ax1.plot(list_t[91:end], list_g[91:end], color="k",
             mfc="none", marker="o", linestyle="");

    # vs T corrected
    ax2.plot(list_t[ 1: 18], list_f[ 1: 18], color="r",
             mfc="none", marker="o", linestyle="");
    ax2.plot(list_t[19: 35], list_f[19: 35], color="y",
             mfc="none", marker="o", linestyle="");
    ax2.plot(list_t[36: 58], list_f[36: 58], color="g",
             mfc="none", marker="o", linestyle="");
    ax2.plot(list_t[59: 74], list_f[59: 74], color="b",
             mfc="none", marker="o", linestyle="");
    ax2.plot(list_t[75: 90], list_f[75: 90], color="c",
             mfc="none", marker="o", linestyle="");
    ax2.plot(list_t[91:end], list_f[91:end], color="k",
             mfc="none", marker="o", linestyle="");

    # genrate theoretical responses
    initialize_vars();
    mod_ts = collect(FT, 19:1:35);
    mod_gs = similar(mod_ts);
    mod_gr = similar(mod_ts);
    for i in eachindex(mod_ts)
        spac.envir.t_air = mod_ts[i] + 273.15;
        spac.envir.p_sat = saturation_vapor_pressure(spac.envir.t_air);
        spac.envir.p_H₂O = spac.envir.p_sat - spac.envir.vpd;
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, FT(0.15));
        mod_gs[i] = gn;
        mod_gr[i] = gn / relative_diffusive_coefficient(can.ps.T);
    end
    ax1.plot(mod_ts, mod_gs, color="silver", linewidth=10, zorder=1);
    ax2.plot(mod_ts, mod_gr, color="silver", linewidth=10, zorder=1);

    # set axes ticks and labels
    set_ylims!(axs, [[0.0,0.12] for i in 1:2]);
    set_xlabels!(axs, [LS_Tleaf_unit for i in 1:2]);
    set_ylabels!(axs, [LS_gwn_unit, LS_gwn25_unit]);
    set_titles!(axs, loc="left");

    # save figure
    fig.set_tight_layout(true);
    if saving
        fig.savefig("figures/2020_nocturnal_gs/4_model_comparison.pdf",
                    bbox_inches="tight");
    end

    #
    #
    # For temperature with higher daytime temperature
    #
    #
    # read the data
    file_t = joinpath(artifact"2020_nocturnal_gs", "g_tem.csv");
    data_t = DataFrame(CSV.File(file_t));
    list_t = (data_t).Tl;
    list_g = (data_t).E ./ (data_t).D .* 86.0;
    list_f = list_g ./ (data_t).F;

    # plot the figure
    fig,axs = create_canvas("GN-S3", ncol=2);
    ax1,ax2 = axs;

    # vs T not corrected
    ax1.plot(list_t[ 1: 18], list_g[ 1: 18], color="r",
             mfc="none", marker="o", linestyle="");
    ax1.plot(list_t[19: 35], list_g[19: 35], color="y",
             mfc="none", marker="o", linestyle="");
    ax1.plot(list_t[36: 58], list_g[36: 58], color="g",
             mfc="none", marker="o", linestyle="");
    ax1.plot(list_t[59: 74], list_g[59: 74], color="b",
             mfc="none", marker="o", linestyle="");
    ax1.plot(list_t[75: 90], list_g[75: 90], color="c",
             mfc="none", marker="o", linestyle="");
    ax1.plot(list_t[91:end], list_g[91:end], color="k",
             mfc="none", marker="o", linestyle="");

    # vs T corrected
    ax2.plot(list_t[ 1: 18], list_f[ 1: 18], color="r",
             mfc="none", marker="o", linestyle="");
    ax2.plot(list_t[19: 35], list_f[19: 35], color="y",
             mfc="none", marker="o", linestyle="");
    ax2.plot(list_t[36: 58], list_f[36: 58], color="g",
             mfc="none", marker="o", linestyle="");
    ax2.plot(list_t[59: 74], list_f[59: 74], color="b",
             mfc="none", marker="o", linestyle="");
    ax2.plot(list_t[75: 90], list_f[75: 90], color="c",
             mfc="none", marker="o", linestyle="");
    ax2.plot(list_t[91:end], list_f[91:end], color="k",
             mfc="none", marker="o", linestyle="");

    # genrate theoretical responses
    initialize_vars();
    mod_ts = collect(FT, 19:1:35);
    mod_gs = similar(mod_ts);
    mod_gr = similar(mod_ts);
    for i in eachindex(mod_ts)
        spac.envir.t_air = mod_ts[i] + 273.15;
        spac.envir.p_sat = saturation_vapor_pressure(spac.envir.t_air);
        spac.envir.p_H₂O = spac.envir.p_sat - spac.envir.vpd;
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, FT(0.15); ΔT=FT(10));
        mod_gs[i] = gn;
        mod_gr[i] = gn / relative_diffusive_coefficient(can.ps.T);
    end
    ax1.plot(mod_ts, mod_gs, color="silver", linewidth=10, zorder=1);
    ax2.plot(mod_ts, mod_gr, color="silver", linewidth=10, zorder=1);

    # set axes ticks and labels
    set_ylims!(axs, [[0.0,0.12] for i in 1:2]);
    set_xlabels!(axs, [LS_Tleaf_unit for i in 1:2]);
    set_ylabels!(axs, [LS_gwn_unit, LS_gwn25_unit]);
    set_titles!(axs, loc="left");

    # save figure
    fig.set_tight_layout(true);
    if saving
        fig.savefig("figures/2020_nocturnal_gs/s3_model_comparison.pdf",
                    bbox_inches="tight");
    end

    #
    #
    # response to respiration rate
    #
    #
    # read data
    file_r = joinpath(artifact"2020_nocturnal_gs", "g_r.csv");
    data   = DataFrame(CSV.File(file_r));

    # create figure and axes
    fig,axs  = create_canvas("GN-5", figsize=(3.7,3.5));
    ax1      = axs[1];
    list_col = ["red", "yellow", "green", "blue", "cyan", "black"];

    # plot data
    for numb in 1:6
        mask = ((data).N .== numb);
        ax1.plot(-(data).A[mask],
                  (data).gswn[mask],
                  marker="o",
                  color=list_col[numb],
                  mfc="none",
                  linestyle="",
                  alpha=0.7);
        plot_line_regress(ax1,
                 -(data).A[mask],
                  (data).gswn[mask],
                  interval=true,
                  color=list_col[numb],
                  alpha=0.3);
    end

    # curve fit to all data
    plot_line_regress(ax1, -(data).A, (data).gswn, interval=true, color="m",
                      alpha=0.3);

    # genrate theoretical responses
    initialize_vars();
    mod_rs = collect(FT, 0.1:0.1:2.0);
    mod_gs = similar(mod_rs);
    for i in eachindex(mod_rs)
        can.ps.Rd25 = mod_rs[i];
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, FT(0.15));
        mod_gs[i] = gn;
    end
    ax1.plot(mod_rs, mod_gs, color="silver", linewidth=10, zorder=1);

    # set labels and limits
    ax1.set_xlim( 0.20,2.0);
    ax1.set_ylim(-0.01,0.4);
    ax1.set_ylabel(LS_gwn_unit, fontsize=16);
    ax1.set_xlabel(LS_Rn_unit  , fontsize=16);

    # save figure
    fig.set_tight_layout(true);
    if saving
        fig.savefig("figures/2020_nocturnal_gs/5_model_comparison.pdf",
                    bbox_inches="tight");
    end

    #
    #
    # response to respiration rate with ΔT
    #
    #
    # create figure and axes
    fig,axs  = create_canvas("GN-S4", figsize=(3.7,3.5));
    ax1      = axs[1];
    list_col = ["red", "yellow", "green", "blue", "cyan", "black"];

    # plot data
    for numb in 1:6
        mask = ((data).N .== numb);
        ax1.plot(-(data).A[mask],
                  (data).gswn[mask],
                  marker="o",
                  color=list_col[numb],
                  mfc="none",
                  linestyle="",
                  alpha=0.7);
        plot_line_regress(ax1,
                 -(data).A[mask],
                  (data).gswn[mask],
                  interval=true,
                  color=list_col[numb],
                  alpha=0.3);
    end

    # curve fit to all data
    plot_line_regress(ax1, -(data).A, (data).gswn, interval=true, color="m",
                      alpha=0.3);

    # genrate theoretical responses
    initialize_vars();
    mod_rs = collect(FT, 0.1:0.1:2.0);
    mod_gs = similar(mod_rs);
    for i in eachindex(mod_rs)
        can.ps.Rd25 = mod_rs[i];
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, FT(0.15); ΔT=FT(10));
        mod_gs[i] = gn;
    end
    ax1.plot(mod_rs, mod_gs, color="silver", linewidth=10, zorder=1);

    # set labels and limits
    ax1.set_xlim( 0.20,2.0);
    ax1.set_ylim(-0.01,0.4);
    ax1.set_ylabel(LS_gwn_unit, fontsize=16);
    ax1.set_xlabel(LS_Rn_unit  , fontsize=16);

    # save figure
    fig.set_tight_layout(true);
    if saving
        fig.savefig("figures/2020_nocturnal_gs/s4_model_comparison.pdf",
                    bbox_inches="tight");
    end

    return nothing
end








###############################################################################
#
# Figure 6: gswn responses to time
#
###############################################################################
"""
    plot_gswn_vs_time(
                project::NocturnalGS2020{FT};
                saving::Bool = false,
                use_latex::Bool=true
    ) where {FT<:AbstractFloat}

Plot gswn vs time, given
- `project` NocturnalGS2020 type project control
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function plot_gswn_vs_time(
            project::NocturnalGS2020{FT};
            saving::Bool = false,
            use_latex::Bool=true
) where {FT<:AbstractFloat}
    # use serif and latex
    if use_latex use_serif_tex(); end

    # read data
    file_t = joinpath(artifact"2020_nocturnal_gs", "g_time.csv");
    data   = DataFrame(CSV.File(file_t));

    # create figure and axes
    fig,axs  = create_canvas("GN-6"; ncol=2, figsize=(7.5,3.7));
    ax1,ax2  = axs;
    tx2      = ax2.twinx();
    list_col = ["red", "yellow", "green", "blue", "cyan", "black"];

    # plot data
    for numb in 1:6
        mask = ((data).N .== numb);
        ax1.plot(-(data).A[mask],
                  (data).gswn[mask],
                  marker="o",
                  color=list_col[numb],
                  mfc="none",
                  linestyle="",
                  alpha=0.3);
        plot_line_regress(ax1,
                 -(data).A[mask],
                  (data).gswn[mask],
                  interval=true,
                  color=list_col[numb],
                  alpha=0.3)
        if numb==1
            lx = collect(1:length((data).gswn[mask])) ./ 2;
            ax2.plot(lx, (data).gswn[mask], "ro", mfc="none", alpha=0.5);
            tx2.plot(lx, -(data).A[mask] , "k-", linewidth=1, alpha=0.5);
        end
    end

    # set labels and limits
    ax1.set_xlim(0.40,1.6);
    ax1.set_ylim(0.02,0.1);
    ax1.set_ylabel(LS_gwn_unit, fontsize=16);
    ax1.set_xlabel(LS_Rn_unit  , fontsize=16);
    ax2.set_xlabel("Time (min)", fontsize=16);
    tx2.set_ylabel(LS_Rn_unit  , fontsize=16);
    set_titles!(axs, labels=String[], loc="left");

    # save figure
    fig.set_tight_layout(true);
    if saving
        fig.savefig("figures/2020_nocturnal_gs/6_g_vs_time.pdf",
                    bbox_inches="tight");
    end

    return nothing
end








###############################################################################
#
# Figure 7: model prediction
#
###############################################################################
"""
    plot_model_extension(
                proj::NocturnalGS2020{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot model extension (VPD and fitness factor change at the same time), given
- `proj` NocturnalGS2020 type project control
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and san-serif font
"""
function plot_model_extension(
            proj::NocturnalGS2020{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    # use sans and latex
    if use_latex use_serif_tex(); end

    # define variables to work with
    spac = create_spac(proj);
    can  = CanopyLayer{FT}(n_leaf=1);
    psm  = C3CLM(FT);
    dsm  = OSMWang{FT}();
    psm.ReT = RespirationTDBernacchi(FT);

    @inline initialize_vars() = (
        spac = create_spac(proj);
        can  = CanopyLayer{FT}(n_leaf=1);
        can.ps.APAR     = 1000;
        spac.envir.wind = 0.1;
    );

    # generate arrays of gain and risk
    list_e = collect(FT, 1:50:2001) * FT(1e-6);
    list_c = similar(list_e);
    list_b = similar(list_e);
    @inline generate_curves() = (
        for i in eachindex(list_e)
            list_b[i],list_c[i] = marginal_gain_risk!(proj, spac, can, psm,
                                                      dsm, list_e[i]);
        end
    );

    # plot the figure
    fig,axs = create_canvas("GN-7"; ncol=2);
    ax1,ax2 = axs;

    # read data
    file_d = joinpath(artifact"2020_nocturnal_gs", "g_vpd.csv");
    file_t = joinpath(artifact"2020_nocturnal_gs", "g_tem.csv");
    data_d = DataFrame(CSV.File(file_d));
    data_t = DataFrame(CSV.File(file_t));
    list_t = (data_t).Tl;
    list_g = (data_t).E ./ (data_t).D .* 86.0;
    list_f = list_g ./ (data_t).F;

    # vs D
    ax1.plot(data_d.D, data_d.Gn, "ko", mfc="none");
    plot_line_regress(ax1, data_d.D, data_d.Gn; interval=true);

    # genrate theoretical responses
    initialize_vars();
    mod_ds = collect(FT, 500:50:3000);
    mod_ff = FT(0.18) .- mod_ds ./ 1000 .* FT(0.03);
    mod_gs = similar(mod_ds);
    for i in eachindex(mod_ds)
        spac.envir.vpd   = mod_ds[i];
        spac.envir.p_H₂O = spac.envir.p_sat - spac.envir.vpd;
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, mod_ff[i]);
        mod_gs[i] = gn;
    end
    ax1.plot(mod_ds ./ 1000, mod_gs, color="silver", linewidth=5, zorder=1,
             label=LS_ff * " \$= 0.18 - 0.03 \\cdot D\$");

    # vs T not corrected
    ax2.plot(list_t, list_g, "ko", mfc="none");
    plot_line_regress(ax2, list_t, list_g; interval=true);

    # genrate theoretical responses
    initialize_vars();
    mod_ts = collect(FT, 19:1:35);
    mod_ff = FT.(0.0235 .* exp.(0.074 .* mod_ts));
    mod_gs = similar(mod_ts);
    for i in eachindex(mod_ts)
        spac.envir.t_air   = mod_ts[i] + 273.15;
        spac.envir.p_sat   = saturation_vapor_pressure(spac.envir.t_air);
        spac.envir.p_H₂O   = spac.envir.p_sat - spac.envir.vpd;
        en,gn = optimal_en_gn(proj, spac, can, psm, dsm, mod_ff[i]);
        mod_gs[i] = gn;
    end
    ax2.plot(mod_ts, mod_gs, color="silver", linewidth=5, zorder=1,
             label=LS_ff * " = \$0.0235 \\cdot \\exp(0.074 \\cdot T)\$");

    # add legend for each axis
    for ax in axs
        ax.legend(loc="lower right");
    end

    set_xlims!(axs, [[0.5, 3.0], [18, 36]]);
    set_ylims!(axs, [[0,0.125], [0.022,0.052]]);
    set_xlabels!(axs, ["D (kPa)", LS_Tleaf_unit]);
    set_ylabels!([ax1], [LS_gwn_unit for i in 1:2]);
    set_titles!(axs, loc="left");

    # save figure
    fig.set_tight_layout(true);
    if saving
        fig.savefig("figures/2020_nocturnal_gs/7_model_extension.pdf",
                bbox_inches="tight");
    end

    return nothing
end








###############################################################################
#
# Figure S1: stable leaf temperature
#
###############################################################################
"""
    plot_si_t_leaf(
                proj::NocturnalGS2020{FT};
                saving::Bool = false,
                use_latex::Bool = true
    ) where {FT<:AbstractFloat}

Plot model extension (VPD and fitness factor change at the same time), given
- `proj` NocturnalGS2020 type project control
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and san-serif font
"""
function plot_si_t_leaf(
            proj::NocturnalGS2020{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    # use serif and latex
    if use_latex use_serif_tex(); end

    # load data
    t_array = [25.002, 24.999, 25.006, 25.004, 25.011, 25.015, 25.027, 25.018,
               25.027, 25.018, 24.991, 24.994, 24.980, 25.001, 24.982, 24.991,
               24.993, 25.005, 25.007, 24.999, 24.997, 25.000, 24.991, 24.999,
               25.015, 25.013, 24.995, 24.998, 24.958, 25.003, 24.993, 25.002,
               25.001, 25.004, 24.986, 25.000, 25.007, 25.014, 25.015, 24.974,
               25.008, 24.993, 24.990, 25.002, 24.992, 24.998, 25.001, 24.997,
               25.007, 24.999, 25.004, 25.003, 24.997, 24.993, 25.006, 25.000,
               25.007, 24.996, 24.998, 24.997, 25.002, 24.990, 24.985, 25.009,
               24.995, 24.998, 25.001, 25.002, 24.994, 24.992, 24.991, 25.006,
               24.997, 24.993, 25.001, 25.002, 24.999, 24.988, 25.002, 24.999,
               24.997, 24.997, 24.990, 24.997, 24.999, 25.002];
    d_array = [0.6346, 0.7928, 1.0441, 1.2289, 1.4186, 1.6108, 1.8164, 2.0285,
               2.2655, 2.3965, 2.2007, 2.0124, 1.8208, 1.6380, 1.4432, 1.2551,
               1.0758, 0.5555, 0.7907, 1.0215, 1.2527, 1.4890, 1.7086, 1.9217,
               2.1316, 2.3461, 0.5856, 0.7986, 0.9160, 0.5960, 0.8500, 1.1118,
               1.3667, 1.6273, 1.8811, 2.1419, 2.4086, 2.6653, 2.7606, 0.5947,
               0.5937, 0.5850, 0.8532, 1.1120, 1.3702, 1.6300, 1.8891, 2.1426,
               2.4083, 2.5975, 2.5071, 0.5869, 0.5891, 0.8564, 1.1138, 1.7170,
               1.6328, 1.8886, 2.1439, 2.2791, 2.1979, 2.1336, 0.5765, 0.5980,
               0.8590, 1.1150, 1.3740, 1.6323, 1.8892, 2.0023, 1.9522, 2.6689,
               2.9236, 0.5882, 0.5920, 0.8578, 1.1155, 1.3706, 1.6296, 1.8907,
               2.1472, 2.4048, 2.6607, 2.1990, 0.5887, 0.6046];

    # create figure and axes
    fig,axs = create_canvas("GN-S1");
    ax1     = axs[1];

    ax1.plot(d_array, t_array, "ko");
    set_xylabels!(axs, ["VPD (kPa)"], [LS_Tleaf_unit]);

    # save figure
    fig.set_tight_layout(true);
    if saving
        fig.savefig("figures/2020_nocturnal_gs/s1_t_leaf.pdf",
                bbox_inches="tight");
    end

    return nothing
end
