###############################################################################
#
# Fit the LeafBio to match spectrum data for Renato's project
#
# To use this script, run this command
# using ResearchProjects; FT=Float64; proj=ClumpingFactor2020{FT}();
# fit_leafbio(proj);
#
###############################################################################
"""
    spectrum_diff!(
                proj::ClumpingFactor2020{FT},
                wls::WaveLengths{FT},
                leaf::LeafBios{FT},
                vars::Array{FT,1},
                refs::Array{FT,1},
                sims::Array{FT,1}
    ) where {FT<:AbstractFloat}

Calculate the sum of square of difference between providing spectrum values and
    simulated spectrum values, given
- `proj` [`ClumpingFactor2020`](@ref) type project control
- `wls` WaveLengths type struct from CanopyLayers module
- `leaf` LeafBio type struct from CanopyLayers module
- `vars` Array of variables to fit
- `refs` Array of reference spectrum data
- `sims` Array of simulated spectrum data (updated within function)
"""
function spectrum_diff!(
            proj::ClumpingFactor2020{FT},
            wls::WaveLengths{FT},
            leaf::LeafBios{FT},
            vars::Array{FT,1},
            refs::Array{FT,1},
            sims::Array{FT,1}
) where {FT<:AbstractFloat}
    leaf.N   = vars[1];
    leaf.Cab = vars[2];
    leaf.Car = vars[3];
    leaf.Ant = vars[4];
    leaf.Cs  = vars[5];
    leaf.Cw  = vars[6];
    leaf.Cm  = vars[7];
    leaf.Cx  = vars[8];
    leaf.fqe = vars[9];

    fluspect!(leaf, wls);

    # pre-allocate them to reduce memory allocation
    mask_par = (wls.WL .> 400) .* (wls.WL .<  700);
    mask_nir = (wls.WL .> 700) .* (wls.WL .< 2500);

    sims[1] = nanmean(leaf.ρ_SW[mask_par]);
    sims[2] = nanmean(leaf.ρ_SW[mask_nir]);
    sims[3] = nanmean(leaf.τ_SW[mask_par]);
    sims[4] = nanmean(leaf.τ_SW[mask_nir]);

    return -sum( (sims .- refs) .^ 2 )
end




"""
    fit_leafbio(proj::ClumpingFactor2020{FT},
                refs::Array{FT,1},
                use_latex::Bool) where {FT<:AbstractFloat}

Fit leaf bio parameters to match given spectrum data, given
- `proj` [`ClumpingFactor2020`](@ref) type project control
- `refs` Array of reference spectrum data
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function fit_leafbio(
            proj::ClumpingFactor2020{FT},
            refs::Array{FT,1} = FT[0.0735, 0.3912, 0.0566, 0.4146];
            saving::Bool = false,
            use_latex::Bool=true
) where {FT<:AbstractFloat}
    # use latex and serif
    if use_latex use_serif_tex(); end

    # 1. create variables to work on
    can  = create_canopy_rt(FT, nLayer=20);
    wls  = create_wave_length(FT);
    dims = create_rt_dims(can, wls);
    leaf = create_leaf_bios(FT, dims);
    sims = similar(refs);
    fluspect!(leaf, wls);

    # 2. fit the data, starting with Renato's initial guess
    x_inis = FT[1.6, 30, 5, 2.75,   0.01, 0.005, 0.001,   0, 0.01];
    #x_inis = FT[1.6, 30,  5, 2.75, 0.01, 0.005, 0.00,   0, 0.01];
    x_mins = FT[1.2, 10,  1,    0,    0, 0.000, 0.00,   0, 0.01];
    x_maxs = FT[2.0, 50, 15,   15, 0.20, 0.050, 0.001,   1, 1.00];
    Δ_inis = FT[0.2, 10,  3,    1, 0.10, 0.005, 0.1, 0.1, 0.10];
    Δ_tole = Δ_inis .* 0.0009;

    # Renato
    #x_inis = FT[1.4, 40, 10, 0,   0, 0.009, 0.012,  0,  1.0];
    #x_mins = FT[1.2, 30,  5, 0,   0, 0.000, 0.000,  0,  0.0];
    #x_maxs = FT[1.6, 50, 15, 5,   0, 0.050,   0.5,  1,  1.0];
    #Δ_inis = FT[0.2, 10,  3, 4, 0.1, 0.005,  0.05, 0.1, 0.1];
    #Δ_inis = Δ_inis/1.e0;
    #Δ_tole = Δ_inis .* 9e-5;

    ms    = ReduceStepMethodND{FT}(x_mins = x_mins,
                                   x_maxs = x_maxs,
                                   x_inis = x_inis,
                                   Δ_inis = Δ_inis);
    st    = SolutionToleranceND{FT}(Δ_tole, 50);
    _f(x) = spectrum_diff!(proj, wls, leaf, x, refs, sims);
    sol   = find_peak(_f, ms, st);
    @show sol;
    spectrum_diff!(proj, wls, leaf, sol, refs, sims);
    @show sum( (sims .- refs) .^ 2 );

    fluspect!(leaf, wls);

    _fig,_axes = create_canvas("CF-1", ncol=2);
    _ax1,_ax2  = _axes;

    _ax1.plot(wls.WL, leaf.τ_SW, "r-", label="\$\\uptau_\\text{leaf}\$");
    _ax1.plot(wls.WL, leaf.ρ_SW, "b-", label="\$\\uprho_\\text{leaf}\$");
    _ax1.set_xlim(350,2410);
    _ax1.set_ylim(0,0.51);
    _ax1.legend(loc="lower right");

    _texts = ["\$\\uprho_\\text{leaf,PAR}\$", "\$\\uprho_\\text{leaf,NIR}\$",
              "\$\\uptau_\\text{leaf,PAR}\$", "\$\\uptau_\\text{leaf,NIR}\$"];
    _ax2.plot(refs[1], sims[1], "bo", mfc="none", label=_texts[1]);
    _ax2.plot(refs[2], sims[2], "bo",             label=_texts[2]);
    _ax2.plot(refs[3], sims[3], "ro", mfc="none", label=_texts[3]);
    _ax2.plot(refs[4], sims[4], "ro",             label=_texts[4]);
    _ax2.plot([0,0.5], [0, 0.5], "k:");
    _ax2.set_xlim(0,0.5);
    _ax2.set_ylim(0,0.5);
    _ax2.legend(loc="lower right");

    # set titles and labels
    set_xlabels!(_axes, ["Wave Length (nm)", "RAMI4PILS reference"]);
    set_ylabels!(_axes,
                 ["\$\\uprho_\\text{leaf}\$ or \$\\uptau_\\text{leaf}\$",
                  "Fluspect model"]);
    set_titles!(_axes; loc="left", paren=false);

    # save figure if you need
    _fig.set_tight_layout(true);

    return sol
end
