using CairoMakie: lines!, surface!
using GeoMakie: Colorbar, Relative, coastlines

using JuliaUtilities.EmeraldRegression: linear_regress
using JuliaUtilities.PlotPlants: plot_hexbin, plot_line_regress!


include("process_data.jl")


# function to plot global maps on axis
function plot_map_on_axis!(ax, cax, lons::Vector, lats::Vector, data::Matrix, clabel::String; vminvmax=nothing)
    lines!(ax, coastlines(), color = :black);
    if isnothing(vminvmax)
        _cm = surface!(ax, lons, lats, data; shading = false);
    else
        _cm = surface!(ax, lons, lats, data; shading = false, colorrange=vminvmax);
    end;
    Colorbar(cax, _cm; height = Relative(0.9), label=clabel);

    return nothing
end;


# function to plot on each axis
function plot_comparison!(fig, ax, ref_1x, clm_1x; laimask=false)
    # plot the hexbin
    _ref_mean, _clm_mean = processed_data(ref_1x, clm_1x; laimask=laimask);
    _cm1 = plot_hexbin(ax, _ref_mean[:], _clm_mean[:], logbins=true);
    plot_line_regress!(ax, (_ref_mean[:],1), _clm_mean[:]; interval=true);
    fig.colorbar(_cm1, ax=ax, label="Number of observations (-)");
    _lr1 = linear_regress((_ref_mean[:],1), _clm_mean[:]);
    _str = "y = $(round(_lr1.COEF[1];digits=3))x";
    _str *= _lr1.COEF[2] > 0 ? " + $(round(_lr1.COEF[2];digits=3))" : " - $(-round(_lr1.COEF[2];digits=3))";
    _str *= "\nR² = $(round(_lr1.R²;digits=3))";

    # plot 1:1 line
    _x_min,_x_max = ax.get_xlim();
    _y_min,_y_max = ax.get_ylim();
    _line_max = min(_x_max, _y_max) - max(_x_max-_x_min, _y_max-_y_min)*0.02;
    _line_min = max(_x_min, _y_min) + max(_x_max-_x_min, _y_max-_y_min)*0.02;
    ax.plot([_line_min, _line_max], [_line_min, _line_max], "k:");
    ax.text(_line_max, _line_min, _str; ha = "right", va = "bottom");

    return nothing
end;


# function to plot on each axis
function plot_comparison_3d!(fig, ax, ref_1x, clm_1x)
    # plot the hexbin
    _ref_mean, _clm_mean = processed_data_3d(ref_1x, clm_1x);
    _cm1 = plot_hexbin(ax, _ref_mean[:], _clm_mean[:], logbins=true);
    plot_line_regress!(ax,( _ref_mean[:],1), _clm_mean[:]; interval=true);
    fig.colorbar(_cm1, ax=ax, label="Number of observations (-)");

    # plot 1:1 line
    _x_min,_x_max = ax.get_xlim();
    _y_min,_y_max = ax.get_ylim();
    _line_max = min(_x_max, _y_max) - max(_x_max-_x_min, _y_max-_y_min)*0.02;
    _line_min = max(_x_min, _y_min) + max(_x_max-_x_min, _y_max-_y_min)*0.02;
    ax.plot([_line_min, _line_max], [_line_min, _line_max], "k:");

    return nothing
end;
