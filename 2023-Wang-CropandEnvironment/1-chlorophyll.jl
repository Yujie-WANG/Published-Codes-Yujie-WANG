using ProgressMeter: @showprogress

using Emerald.EmeraldLand.Namespace: MultiLayerSPAC, SPACConfiguration
using Emerald.EmeraldLand.SPAC: CNPP, GPP, PPAR, initialize!, soil_plant_air_continuum!, update!
using Emerald.EmeraldIO.Netcdf: append_nc!, create_nc!, read_nc

using Emerald.EmeraldVisualization: canvas, decorate!, save_canvas!


# mkpath if the path does not exist
if !isdir("$(@__DIR__)/output")
    mkpath("$(@__DIR__)/output");
end;
OUTPUT_FOLDER = "$(@__DIR__)/output";


# function to compute CNPP based on LAI and chlorophyll content
function cnpp(lai::Number, chl::Number; FT::DataType = Float64)
    _config = SPACConfiguration{FT}();
    _spac = MultiLayerSPAC(_config);
    initialize!(_spac, _config);
    update!(_spac, _config; cab = chl, car = chl / 7, lai = lai, vcmax_expo = 0.15);
    for _i in 1:60
        soil_plant_air_continuum!(_spac, _config, FT(60); p_on = false, t_on = false, θ_on = false);
    end;
    _cost = lai * _spac.LEAVES[1].BIO.lma * 2 * 1e4 / 30 * 1e6 / (4 * 30 * 6 * 3600);

    return CNPP(_spac) - _cost
end


# function to generate data to plot
function generate_data!()
    # data for contour
    _lais = collect(Float64, 1:0.25:6);
    _chls = collect(Float64, 2:1:80);
    _matr = zeros(Float64, length(_chls), length(_lais));
    @showprogress for _i in eachindex(_chls), _j in eachindex(_lais)
        _matr[_i,_j] = cnpp(_lais[_j], _chls[_i]);
    end;
    create_nc!("$(OUTPUT_FOLDER)/1-chlorophyll.nc", String["chl", "lai"], [length(_chls), length(_lais)]);
    append_nc!("$(OUTPUT_FOLDER)/1-chlorophyll.nc", "chl", _chls, Dict{String,String}("about" => "Chlorophyll content"), String["chl"]);
    append_nc!("$(OUTPUT_FOLDER)/1-chlorophyll.nc", "lai", _lais, Dict{String,String}("about" => "Leaf area index"), String["lai"]);
    append_nc!("$(OUTPUT_FOLDER)/1-chlorophyll.nc", "npp", _matr, Dict{String,String}("about" => "Canopy NPP"), String["chl", "lai"]);

    return nothing
end


# function to compute GPP and PPAR based on LAI and chlorophyll content
function gpp_ppar(lai::Number, chl::Number; FT::DataType = Float64)
    _config = SPACConfiguration{FT}();
    _spac = MultiLayerSPAC(_config);
    initialize!(_spac, _config);
    update!(_spac, _config; cab = chl, car = chl / 7, lai = lai, vcmax_expo = 0.15);
    for _i in 1:60
        soil_plant_air_continuum!(_spac, _config, FT(60); p_on = false, t_on = false, θ_on = false);
    end;

    return GPP(_spac), PPAR(_spac)
end

# function to display some information related GPP and PPAR
function make_examples!()
    @show gpp_ppar(4, 60);
    @show gpp_ppar(4, 35);
    @show gpp_ppar(4, 10);

    return nothing
end


# function to plot the figure
function plot_optimal_chl!()
    _fig = canvas(1; dpi = 600, figsize = (4.5,3.5));
    _axs = _fig.axes;

    # plot contour
    _chl = read_nc("$(OUTPUT_FOLDER)/1-chlorophyll.nc", "chl");
    _lai = read_nc("$(OUTPUT_FOLDER)/1-chlorophyll.nc", "lai");
    _npp = read_nc("$(OUTPUT_FOLDER)/1-chlorophyll.nc", "npp");
    _opt = [_chl[findmax(_npp[:,_i])[2]] for _i in eachindex(_lai)];
    _cm = _axs[1].contourf(_lai, _chl, _npp);
    _axs[1].plot(_lai, _opt, "k:");
    _fig.colorbar(_cm; ax = _axs[1], label = "Yield (μmol m⁻² s⁻¹)");

    decorate!(_axs; add_title = false, xaxis_labels = "LAI (-)", yaxis_labels = "chl (μg cm⁻²)");
    save_canvas!(_fig, "1-chlorophyll"; folder = OUTPUT_FOLDER, formats = ["pdf","png"]);

    return nothing
end
