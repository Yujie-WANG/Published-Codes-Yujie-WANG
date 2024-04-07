using PkgUtility: ρg_MPa
using PlantHydraulics: LeafHydraulics, LogisticSingle, RootHydraulics, StemHydraulics, TreeSimple, critical_flow, end_pressure, vc_integral, xylem_k_ratio
using PlotPlants: create_canvas, save_canvas!, set_titles!, set_xticks!, set_xylabels!, set_xylims!, use_serif_tex


# Use LaTeX text render. You need PdfLaTeX installed to plot the figures
use_serif_tex();


# Function to plot Figure 2
function plot_measures()
    # define a hydraulic segment, say a leaf hydraulics
    _shs = StemHydraulics{Float64}(N=10, Δh=10);

    # plot the figure
    _fig,_axs = create_canvas(1; figsize=(4,3.5));
    _ax1 = _axs[1];
    _tx1 = _ax1.twinx();

    # plot the supply curve on ax1
    _pdos = collect(Float64, -0.5:-0.1:-3);
    _es   = vc_integral.([_shs.vc], _pdos, -0.5) .* _shs.k_max;
    _ax1.plot(_pdos, _es, "k-", label="\$E\$")

    # plot the vulnerability curve on tx1
    _ps = collect(Float64, 0:-0.1:-3);
    _ks = xylem_k_ratio.([_shs.vc], _ps) .* _shs.k_max;
    _tx1.plot(_ps, _ks, "k--", label="\$k\$");

    # use -2.0 as the base line
    _ax1.plot([-2,-2], [0,75], "k:", alpha=0.5);

    # annotate the quantities
    _slope = (_es[17] - _es[15]) / (_pdos[17] - _pdos[15]);
    _dedp_ps = collect(Float64, -2.3:0.05:-1.6);
    _dedp_ls = _es[16] .- _pdos[16] .* _slope .+ _dedp_ps .* _slope;
    _ax1.plot(_dedp_ps, _dedp_ls, "r:");
    _val_1 = round(_slope/_shs.k_max; digits=2);
    _ax1.annotate("\$\\dfrac{\\text{d}E}{\\text{d}P_\\text{canopy}}\$", xy=(-2,_es[16]), xytext=(-1.5,67), arrowprops=Dict(["arrowstyle" => "->" , "color" => "r"]), color="r");
    _ax1.plot([-0.5,-2.0], [0,_es[16]], "c:");
    _val_2 = round(_es[16]/1.5/_shs.k_max; digits=2);
    _ax1.annotate("\$\\dfrac{E}{\\Delta P}\$", xy=(-1.25,_es[16]/2), xytext=(-1.6,14), arrowprops=Dict(["arrowstyle" => "->" , "color" => "c"]), color="c");
    _tx1.plot(-2, _ks[21], "bo");
    _val_3 = round(_ks[21]/_shs.k_max; digits=2);
    _tx1.annotate("\$k\$", xy=(-2,_ks[21]), xytext=(-2.5,22), arrowprops=Dict(["arrowstyle" => "->" , "color" => "b"]), color="b");

    # add legends
    _ax1.legend(loc="upper left");
    _tx1.legend(loc="upper right");

    # set up the x and y axes
    set_xylabels!([_ax1,_tx1], "\$P\$, \$P_\\text{canopy}\$ (MPa)", ["\$E\$ (mol s\$^{-1}\$)", "\$k\$ (mol MPa\$^{-1}\$ s\$^{-1}\$)"])
    set_xylims!([_ax1,_tx1], [-3,0], [[0,77], [0,58]]);

    # save the figure
    save_canvas!(_fig, "2_quantities.pdf", true);

    return _fig;
end;


# Plot Figure 2 if the figure does not exist
if !isfile("2_quantities.pdf")
    plot_measures();
end


# Function to create a hydraulic system for simulation
function simple_tree()
    # initialize the tree
    _tree = TreeSimple{Float64}(root=RootHydraulics{Float64}(N=10000), stem=StemHydraulics{Float64}(N=10000,Δh=10), leaf=LeafHydraulics{Float64}(N=10000,area=50/0.04));

    # use traits from water birch
    _tree.root.vc.b = 1.879;
    _tree.root.vc.c = 2.396;
    _tree.stem.vc.b = 2.238;
    _tree.stem.vc.c = 9.380;
    _tree.leaf.vc.b = 1.897;
    _tree.leaf.vc.c = 2.203;

    return _tree
end;


# Function to calculate the comparison matrices
function comparison_matrices()
    # initialize the tree
    _tree = simple_tree();

    # create 5 matrices to compare
    _ps = collect(Float64, 0:-0.05:-2.7);
    _es = collect(Float64, 0:0.02:1-eps());
    _mat_dedp = ones(length(_ps),length(_es)) .* NaN;
    _mat_k_r  = ones(length(_ps),length(_es)) .* NaN;
    _mat_k_s  = ones(length(_ps),length(_es)) .* NaN;
    _mat_k_l  = ones(length(_ps),length(_es)) .* NaN;
    _mat_tree = ones(length(_ps),length(_es)) .* NaN;

    # fill the matrices with simulations
    for _i in eachindex(_ps)
        @info "Soil water potential is $(_ps[_i]) MPa.";
        _tree.root.p_ups = _ps[_i];
        _e_crit = critical_flow(_tree);
        _p_zero = end_pressure(_tree, 0.0);
        for _j in eachindex(_es)
            _p0   = end_pressure(_tree, _es[_j] * _e_crit);
            _p1   = end_pressure(_tree, (_es[_j] + 1e-4) * _e_crit);
            _dedp = 1e-4 * _e_crit / (_p0 - _p1);
            _k_r  = xylem_k_ratio(_tree.root.vc, _p0);
            _k_s  = xylem_k_ratio(_tree.stem.vc, _p0);
            _k_l  = xylem_k_ratio(_tree.leaf.vc, _p0);
            _k_t  = (_es[_j] + 1e-4) * _e_crit / (_p_zero - _p1);
            _mat_dedp[_i,_j] = _dedp;
            _mat_k_r[_i,_j]  = _k_r;
            _mat_k_s[_i,_j]  = _k_s;
            _mat_k_l[_i,_j]  = _k_l;
            _mat_tree[_i,_j] = _k_t;
        end;
    end;

    return _ps, _es, _mat_dedp, _mat_k_r, _mat_k_s, _mat_k_l, _mat_tree
end;


# Function to plot Figure 3
function plot_comparison()
    # calculate the matrices
    PS, ES, MAT_DEDP, MAT_K_R, MAT_K_S, MAT_K_L, MAT_TREE = comparison_matrices();

    # rescale tha matrices
    MAT_DEDP ./= maximum(MAT_DEDP);
    MAT_K_R  ./= maximum(MAT_K_R );
    MAT_K_S  ./= maximum(MAT_K_S );
    MAT_K_L  ./= maximum(MAT_K_L );
    MAT_TREE ./= maximum(MAT_TREE);

    # create canvas
    _fig,_axs = create_canvas(2; nrow=2, ncol=5);
    _ax1,_ax2,_ax3,_ax4,_ax5,_bx1,_bx2,_bx3,_bx4,_bx5 = _axs;

    # plot the data
    _cm1 = _ax1.pcolor(ES*100, PS, MAT_DEDP, vmin=0, vmax=1, cmap="coolwarm_r");
    _cm2 = _ax2.pcolor(ES*100, PS, MAT_K_R , vmin=0, vmax=1, cmap="coolwarm_r");
    _cm3 = _ax3.pcolor(ES*100, PS, MAT_K_S , vmin=0, vmax=1, cmap="coolwarm_r");
    _cm4 = _ax4.pcolor(ES*100, PS, MAT_K_L , vmin=0, vmax=1, cmap="coolwarm_r");
    _cm5 = _ax5.pcolor(ES*100, PS, MAT_TREE, vmin=0, vmax=1, cmap="coolwarm_r");

    # rescale tha matrices
    for _i in eachindex(PS)
        MAT_DEDP[_i,:] ./= MAT_DEDP[_i,1];
        MAT_K_R[_i,:]  ./= MAT_K_R[_i,1];
        MAT_K_S[_i,:]  ./= MAT_K_S[_i,1];
        MAT_K_L[_i,:]  ./= MAT_K_L[_i,1];
        MAT_TREE[_i,:] ./= MAT_TREE[_i,1];
    end;

    # plot the data
    _dm1 = _bx1.pcolor(ES*100, PS, MAT_DEDP, vmin=0, vmax=1, cmap="coolwarm_r");
    _dm2 = _bx2.pcolor(ES*100, PS, MAT_K_R , vmin=0, vmax=1, cmap="coolwarm_r");
    _dm3 = _bx3.pcolor(ES*100, PS, MAT_K_S , vmin=0, vmax=1, cmap="coolwarm_r");
    _dm4 = _bx4.pcolor(ES*100, PS, MAT_K_L , vmin=0, vmax=1, cmap="coolwarm_r");
    _dm5 = _bx5.pcolor(ES*100, PS, MAT_TREE, vmin=0, vmax=1, cmap="coolwarm_r");

    # plot the color bar
    _fig.colorbar(_cm1, ax=_ax1);
    _fig.colorbar(_cm2, ax=_ax2);
    _fig.colorbar(_cm3, ax=_ax3);
    _fig.colorbar(_cm4, ax=_ax4);
    _fig.colorbar(_cm5, ax=_ax5);
    _fig.colorbar(_dm1, ax=_bx1);
    _fig.colorbar(_dm2, ax=_bx2);
    _fig.colorbar(_dm3, ax=_bx3);
    _fig.colorbar(_dm4, ax=_bx4);
    _fig.colorbar(_dm5, ax=_bx5);

    # set up the axis
    _s_dedp = "\\dfrac{\\text{d}E}{\\text{d}P_\\text{canopy}}";
    _s_k_km = "\$k : k_\\text{max}\$";
    _s_k_ks = "\$k : k(\\Psi_\\text{soil})\$";
    _labels = ["\$$(_s_dedp) : \\text{max}($(_s_dedp))\$", "Root $(_s_k_km)", "Stem $(_s_k_km)", "Leaf $(_s_k_km)", "Tree $(_s_k_km)",
               "\$$(_s_dedp) : ($(_s_dedp))_{E = 0}\$", "Root $(_s_k_ks)", "Stem $(_s_k_ks)", "Leaf $(_s_k_ks)", "Tree \$k : k_{E = 0}\$"];
    set_xylabels!(_axs, repeat(["","\$E / E_\\text{crit}\$ (\\%)"]; inner=5), repeat(["\$\\Psi_\\text{soil}\$ (MPa)","","","",""]; outer=2));
    set_titles!(_axs; loc="left", labels=_labels);

    # save the figure
    save_canvas!(_fig, "3_comparison.pdf", true);

    return _fig
end;


# Plot Figure 3 if the figure does not exist
if !isfile("3_comparison.pdf")
    plot_comparison();
end


# Function to plot Figure 4
function plot_analytic_estimation()
    # define the stems with height change or no height change
    _stem_0 = StemHydraulics{Float64}(vc=LogisticSingle{Float64}(64,3), k_max = 20, N = 10000, Δh=0);
    _stem_1 = StemHydraulics{Float64}(vc=LogisticSingle{Float64}(64,3), k_max = 20, N = 10000, Δh=5);
    _stem_2 = StemHydraulics{Float64}(vc=LogisticSingle{Float64}(64,3), k_max = 20, N = 10000, Δh=10);
    _stem_3 = StemHydraulics{Float64}(vc=LogisticSingle{Float64}(64,3), k_max = 20, N = 10000, Δh=20);

    # plot the figures
    _fig,_axs = create_canvas(3; ncol=4);
    _ax1,_ax2,_ax3,_ax4 = _axs;

    # plot the case of 0 MPa
    _p_upss = [0.0, -0.5, -1.0, -1.5];
    for _i in eachindex(_p_upss)
        _p_ups = _p_upss[_i];
        _stem_0.p_ups = _p_ups;
        _stem_1.p_ups = _p_ups;
        _stem_2.p_ups = _p_ups;
        _stem_3.p_ups = _p_ups;
        _ps_0  = end_pressure(_stem_0, 0.0);
        _ps_1  = end_pressure(_stem_1, 0.0);
        _ps_2  = end_pressure(_stem_2, 0.0);
        _ps_3  = end_pressure(_stem_3, 0.0);
        _pss_0 = collect(_ps_0:-0.01:-3.01);
        _pss_1 = collect(_ps_1:-0.01:-3.01);
        _pss_2 = collect(_ps_2:-0.01:-3.01);
        _pss_3 = collect(_ps_3:-0.01:-3.01);
        _est_0 = _stem_0.k_max .* vc_integral.([_stem_0.vc], _pss_0, _p_ups) .* (_pss_0 .+ ρg_MPa() .* _stem_0.Δh .- _p_ups) ./ (_pss_0 .- _p_ups);
        _est_1 = _stem_1.k_max .* vc_integral.([_stem_0.vc], _pss_1, _p_ups) .* (_pss_1 .+ ρg_MPa() .* _stem_1.Δh .- _p_ups) ./ (_pss_1 .- _p_ups);
        _est_2 = _stem_2.k_max .* vc_integral.([_stem_0.vc], _pss_2, _p_ups) .* (_pss_2 .+ ρg_MPa() .* _stem_2.Δh .- _p_ups) ./ (_pss_2 .- _p_ups);
        _est_3 = _stem_3.k_max .* vc_integral.([_stem_0.vc], _pss_3, _p_ups) .* (_pss_3 .+ ρg_MPa() .* _stem_3.Δh .- _p_ups) ./ (_pss_3 .- _p_ups);
        _es_0  = vc_integral.([_stem_0.vc], _pss_0, _p_ups, _stem_0.Δh, _est_0, _stem_0.k_max);
        _es_1  = vc_integral.([_stem_1.vc], _pss_1, _p_ups, _stem_1.Δh, _est_1, _stem_1.k_max);
        _es_2  = vc_integral.([_stem_2.vc], _pss_2, _p_ups, _stem_2.Δh, _est_2, _stem_2.k_max);
        _es_3  = vc_integral.([_stem_3.vc], _pss_3, _p_ups, _stem_3.Δh, _est_3, _stem_3.k_max);
        _es_0[1] = 0;
        _es_1[1] = 0;
        _es_2[1] = 0;
        _es_3[1] = 0;
        _axs[_i].plot(_pss_0, _es_0, "k-" , label="\$\\Delta h\$ = 0 m");
        _axs[_i].plot(_pss_1, _es_1, "k--", label="\$\\Delta h\$ = 5 m");
        _axs[_i].plot(_pss_2, _es_2, "k-.", label="\$\\Delta h\$ = 10 m");
        _axs[_i].plot(_pss_3, _es_3, "k:" , label="\$\\Delta h\$ = 20 m");
    end;
    _ax1.legend(loc="lower left");

    # set up the axis
    set_xylabels!(_axs, "\$P_\\text{canopy}\$ (MPa)", ["\$E\$ (mol s\$^{-1}\$)","","",""]);
    set_xylims!(_axs, [[-3,0], [-3,-0.5], [-3,-1], [-3,-1.5]], [[0,30], [0,20], [0,10], [0,4]]);
    set_xticks!([_ax1], [-3,-2.5,-2,-1.5,-1,-0.5,0]);
    set_titles!(_axs; loc="left", labels=["\$\\Psi_\\text{soil} = 0\$ MPa", "\$\\Psi_\\text{soil} = -0.5\$ MPa", "\$\\Psi_\\text{soil} = -1\$ MPa", "\$\\Psi_\\text{soil} = -1.5\$ MPa"]);

    # save the figure
    save_canvas!(_fig, "4_height.pdf", true);

    return _fig
end;


# Plot Figure 4 if the figure does not exist
if !isfile("4_height.pdf")
    close("all");
    plot_analytic_estimation();
end


# Function to plot Figure 5
function plot_surface()
    # define the stems with height change or no height change
    _stem = StemHydraulics{Float64}(vc=LogisticSingle{Float64}(64,3), k_max = 20, N = 10000, Δh=20);

    # create a matrix to save the k results
    _p_upss = collect(Float64, 0:-0.01:-2.5);
    _p_doss = collect(Float64, 0:-0.01:-3.5);
    _mat_k  = ones(length(_p_upss), length(_p_doss)) .* NaN;

    # update the matrix
    for _i in eachindex(_p_upss)
        _p_ups = _p_upss[_i];
        _stem.p_ups = _p_ups;
        _p_0 = end_pressure(_stem, 0.0);
        for _j in eachindex(_p_doss)
            _p_dos = _p_doss[_j];
            if _p_dos < _p_0
                _est = _stem.k_max * vc_integral(_stem.vc, _p_dos, _p_ups) * (_p_dos + ρg_MPa() * _stem.Δh - _p_ups) / (_p_dos - _p_ups);
                _e   = vc_integral(_stem.vc, _p_dos, _p_ups, _stem.Δh, _est, _stem.k_max);
                _k   = _e / (_p_0 - _p_dos);
            else
                _est = _stem.k_max * vc_integral(_stem.vc, _p_0 - 1e-6, _p_ups) * (_p_0 - 1e-6 + ρg_MPa() * _stem.Δh - _p_ups) / (_p_0 - 1e-6 - _p_ups);
                _e   = vc_integral(_stem.vc, _p_0 - 1e-6, _p_ups, _stem.Δh, _est, _stem.k_max);
                _k   = _e * 1e6;
            end;
            _mat_k[_i,_j] = _k;
        end;
    end;

    # create the canvas and axis
    _fig,_axs = create_canvas(4; figsize=(4.4,3.5));
    _ax1, = _axs;

    # plot the figure
    _cm1 = _ax1.contourf(_p_upss, _p_doss, _mat_k');
    _fig.colorbar(_cm1, ax=_ax1, label="\$k_\\text{plant}\$ (mol s\$^{-1}\$ MPa\$^{-1}\$)");
    _ax1.plot([-1.5,0], [-1.7,-3.2], "k:");

    # set up the axis
    set_xylabels!(_axs, "\$\\Psi_\\text{soil}\$ (MPa)", "\$P_\\text{canopy}\$ (MPa)");

    # save the figure
    save_canvas!(_fig, "5_surface.pdf", true);

    return _fig
end;


# Plot Figure 5 if the figure does not exist
if !isfile("5_surface.pdf")
    plot_surface();
end
