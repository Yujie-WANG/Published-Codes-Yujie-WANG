using ProgressMeter: @showprogress

using Emerald.EmeraldLand.LeafOptics: leaf_spectra!
using Emerald.EmeraldLand.Namespace: MultiLayerSPAC, SPACConfiguration
using Emerald.EmeraldLand.SPAC: CNPP, GPP, PPAR, initialize!, soil_plant_air_continuum!, update!
using Emerald.EmeraldIO.Netcdf: append_nc!, create_nc!, read_nc

using Emerald.EmeraldVisualization: canvas, decorate!, save_canvas!


# mkpath if the path does not exist
if !isdir("$(@__DIR__)/output")
    mkpath("$(@__DIR__)/output");
end;
OUTPUT_FOLDER = "$(@__DIR__)/output";


# function to compute sum PPAR and A based on LAI and ci
function gpp_ppar_2(lai::Number, ci::Number; FT::DataType = Float64)
    _config = SPACConfiguration{FT}();
    _spac = MultiLayerSPAC(_config);
    initialize!(_spac, _config);
    update!(_spac, _config; lai = lai, vcmax_expo = 0.15);
    _spac.CANOPY.ci = ci;
    _spac.CANOPY.Ω_A = ci;
    for _ in 1:60
        soil_plant_air_continuum!(_spac, _config, FT(60); p_on = false, t_on = false, θ_on = false);
    end;

    return GPP(_spac), PPAR(_spac)
end


# function to give the number of a few scenarios
function test_cases!()
    @show gpp_ppar_2(3, 1);
    @show gpp_ppar_2(4, 1);
    @show gpp_ppar_2(3, 0.7);

    return nothing
end


# function to compute CNPP based on LAI and ci
function cnpp_2(lai::Number, ci::Number; FT::DataType = Float64)
    _config = SPACConfiguration{FT}();
    _spac = MultiLayerSPAC(_config);
    initialize!(_spac, _config);
    update!(_spac, _config; lai = lai, vcmax_expo = 0.15);
    _spac.CANOPY.ci = ci;
    _spac.CANOPY.Ω_A = ci;
    for _ in 1:60
        soil_plant_air_continuum!(_spac, _config, FT(60); p_on = false, t_on = false, θ_on = false);
    end;
    _cost = lai * _spac.LEAVES[1].BIO.lma * 2 * 1e4 / 30 * 1e6 / (4 * 30 * 6 * 3600);

    return CNPP(_spac) - _cost
end


# function to generate data to plot
function generate_data_2!()
    # data for contour
    _lais = collect(Float64, 1:0.25:6);
    _clis = collect(Float64, 0.4:0.02:1.0);
    _matr = zeros(Float64, length(_clis), length(_lais));
    @showprogress for _i in eachindex(_clis), _j in eachindex(_lais)
        _matr[_i,_j] = cnpp_2(_lais[_j], _clis[_i]);
    end;
    create_nc!("$(OUTPUT_FOLDER)/3-lai-ci.nc", String["cli", "lai"], [length(_clis), length(_lais)]);
    append_nc!("$(OUTPUT_FOLDER)/3-lai-ci.nc", "cli", _clis, Dict{String,String}("about" => "Clumping index"), String["cli"]);
    append_nc!("$(OUTPUT_FOLDER)/3-lai-ci.nc", "lai", _lais, Dict{String,String}("about" => "Leaf area index"), String["lai"]);
    append_nc!("$(OUTPUT_FOLDER)/3-lai-ci.nc", "npp", _matr, Dict{String,String}("about" => "Canopy NPP"), String["cli", "lai"]);

    return nothing
end


# function to plot the figure
function plot_optimal_ci!()
    _fig = canvas(3; dpi = 600, figsize = (4.5,3.5));
    _axs = _fig.axes;

    # plot contour
    _cli = read_nc("$(OUTPUT_FOLDER)/3-lai-ci.nc", "cli");
    _lai = read_nc("$(OUTPUT_FOLDER)/3-lai-ci.nc", "lai");
    _npp = read_nc("$(OUTPUT_FOLDER)/3-lai-ci.nc", "npp");
    _opt = [_cli[findmax(_npp[:,_i])[2]] for _i in eachindex(_lai)];
    _cm = _axs[1].contourf(_lai, _cli, _npp);
    _axs[1].plot(_lai, _opt, "k:");
    _fig.colorbar(_cm; ax = _axs[1], label = "Yield (μmol m⁻² s⁻¹)");

    decorate!(_axs; add_title = false, xaxis_labels = "LAI (-)", yaxis_labels = "CI (-)");
    save_canvas!(_fig, "3-lai-ci"; folder = OUTPUT_FOLDER, formats = ["pdf","png"]);

    return nothing
end


# function to test the multi-farming
function test_multi_farming!()
    # default setup
    _config = SPACConfiguration{Float64}();
    _spac = MultiLayerSPAC(_config);
    initialize!(_spac, _config);
    update!(_spac, _config; lai = 3, vcmax_expo = 0.15);
    for _ in 1:60
        soil_plant_air_continuum!(_spac, _config, 60.0; p_on = false, t_on = false, θ_on = false);
    end;
    @info "Default" GPP(_spac) PPAR(_spac);

    # less greener upper canopy
    _config = SPACConfiguration{Float64}();
    _spac = MultiLayerSPAC(_config);
    for _i in 1:10
        _spac.LEAVES[_i].BIO.cab = 10;
    end;
    initialize!(_spac, _config);
    update!(_spac, _config; lai = 3, vcmax_expo = 0.15);
    for _ in 1:60
        soil_plant_air_continuum!(_spac, _config, 60.0; p_on = false, t_on = false, θ_on = false);
    end;
    @info "Less greener upper canopy" GPP(_spac) PPAR(_spac);

    # lower canopy with NIR mutant
    _config = SPACConfiguration{Float64}();
    _spac = MultiLayerSPAC(_config);
    initialize!(_spac, _config);
    update!(_spac, _config; lai = 3, vcmax_expo = 0.15);
    # change the shape of the absorption feature of chl
    _wls = _config.WLSET;
    _lha = deepcopy(_config.LHA);
    for _i in eachindex(_wls.Λ)[end:-1:1]
        if 680 < _wls.Λ[_i] < 775
            _lha.K_CAB[_i] = _lha.K_CAB[_i-5];
        end;
    end;
    for _i in 7:12
        _spac.LEAVES[_i].BIO._v_storage = 0;
        leaf_spectra!(_spac.LEAVES[_i].BIO, _wls, _lha, _spac.LEAVES[_i].HS.v_storage);
    end;
    for _ in 1:60
        soil_plant_air_continuum!(_spac, _config, 60.0; p_on = false, t_on = false, θ_on = false);
    end;
    @info "NIR mutant lower canopy" GPP(_spac) PPAR(_spac);

    return nothing
end
