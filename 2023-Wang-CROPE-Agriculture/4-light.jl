using Emerald.EmeraldLand.Namespace: HyperspectralRadiation, MultiLayerSPAC, SPACConfiguration
using Emerald.EmeraldLand.SPAC: CNPP, initialize!, soil_plant_air_continuum!, update!

using Emerald.EmeraldVisualization: canvas, decorate!, save_canvas!


# mkpath if the path does not exist
if !isdir("$(@__DIR__)/output")
    mkpath("$(@__DIR__)/output");
end;
OUTPUT_FOLDER = "$(@__DIR__)/output";


# function to get CNPP from lai and light source (on a 12 hour basis)
function cnpp_4(lai::Number, chl::Number, srad::HyperspectralRadiation{FT}) where {FT<:AbstractFloat}
    _config = SPACConfiguration{FT}();
    _spac = MultiLayerSPAC(_config);
    _spac.METEO.rad_sw = deepcopy(srad);
    initialize!(_spac, _config);
    update!(_spac, _config; cab = chl, car = chl / 7, lai = lai, vcmax_expo = 0.15);
    for _ in 1:10
        soil_plant_air_continuum!(_spac, _config, FT(360); p_on = false, t_on = false, θ_on = false);
    end;
    soil_plant_air_continuum!(_spac, _config, FT(0); p_on = false, t_on = false, θ_on = false);
    _cost = lai * _spac.LEAVES[1].BIO.lma * 2 * 1e4 / 30 * 1e6 / (12 * 30 * 6 * 3600);

    return CNPP(_spac) - _cost
end


# function to get the curve of CNPP ~ Energy
function cnpp_light_curve(light_source::String = "nature", lai::Number = 3, chl::Number = 50)
    @info "Working on $(light_source) light with LAI = $(lai) and CHL = $(chl)...";

    _config = SPACConfiguration{Float64}();
    _spac = MultiLayerSPAC(_config);
    _srad_bak = deepcopy(_spac.METEO.rad_sw);
    _wls = _config.WLSET;
    _srad_nature = deepcopy(_srad_bak);
    _srad_green = deepcopy(_srad_bak);
    _srad_red = deepcopy(_srad_bak);
    _srad_nir = deepcopy(_srad_bak);
    _srad_green.e_direct .= 0;
    _srad_green.e_diffuse[.!(500 .< _wls.Λ .< 570)] .= 0;
    _srad_red.e_direct .= 0;
    _srad_red.e_diffuse[.!(620 .< _wls.Λ .< 700)] .= 0;
    _srad_nir.e_direct .= 0;
    _srad_nir.e_diffuse[.!(620 .< _wls.Λ .< 730)] .= 0;

    _es = Float64[];
    _ps = Float64[];
    for _rad in 25:25:600
        _srad = if light_source == "nature"
            deepcopy(_srad_nature)
        elseif light_source == "green"
            deepcopy(_srad_green)
        elseif light_source == "red"
            deepcopy(_srad_red)
        elseif light_source == "nir"
            deepcopy(_srad_nir)
        else
            error("Light source $(light_source) not supported!")
        end;
        _srad_sum = _wls.ΔΛ' * (_srad.e_direct .+ _srad.e_diffuse) * 1e-3
        _ratio = _rad / _srad_sum;
        _srad.e_direct .*= _ratio;
        _srad.e_diffuse .*= _ratio;
        push!(_es, _rad);
        push!(_ps, cnpp_4(lai, chl, _srad));
    end;

    return _es, _ps
end


# function to get the curve of CNPP ~ Energy
function optimal_light_intensity(light_source::String = "nature", lai::Number = 3, chl::Number = 50)
    @info "Working on $(light_source) light with LAI = $(lai) and CHL = $(chl)...";

    _config = SPACConfiguration{Float64}();
    _spac = MultiLayerSPAC(_config);
    _srad_bak = deepcopy(_spac.METEO.rad_sw);
    _wls = _config.WLSET;
    _srad_nature = deepcopy(_srad_bak);
    _srad_green = deepcopy(_srad_bak);
    _srad_red = deepcopy(_srad_bak);
    _srad_nir = deepcopy(_srad_bak);
    _srad_green.e_direct .= 0;
    _srad_green.e_diffuse[.!(500 .< _wls.Λ .< 570)] .= 0;
    _srad_red.e_direct .= 0;
    _srad_red.e_diffuse[.!(620 .< _wls.Λ .< 700)] .= 0;
    _srad_nir.e_direct .= 0;
    _srad_nir.e_diffuse[.!(620 .< _wls.Λ .< 730)] .= 0;

    _max = -Inf;
    _max_rad = 0;
    for _rad in 10:10:600
        _srad = if light_source == "nature"
            deepcopy(_srad_nature)
        elseif light_source == "green"
            deepcopy(_srad_green)
        elseif light_source == "red"
            deepcopy(_srad_red)
        elseif light_source == "nir"
            deepcopy(_srad_nir)
        else
            error("Light source $(light_source) not supported!")
        end;
        _srad_sum = _wls.ΔΛ' * (_srad.e_direct .+ _srad.e_diffuse) * 1e-3
        _ratio = _rad / _srad_sum;
        _srad.e_direct .*= _ratio;
        _srad.e_diffuse .*= _ratio;
        _yield = cnpp_4(lai, chl, _srad);
        if _yield / _rad > _max
            _max = _yield / _rad;
            _max_rad = _rad
        else
            break;
        end;
    end;

    return _max,_max_rad
end


# function to plot the optimal light source
function plot_light_source!()
    # create a canvas to plot the figure
    _fig = canvas("4"; dpi = 600, ncol = 3, nrow = 2, figsize = (12,6.5));
    _axs = _fig.axes;

    # plot nature light on axis 1
    _es,_ps = cnpp_light_curve("nature", 2);
    _max = maximum(_ps ./ _es);
    _axs[1].plot(_es, _ps, "k:"; label = "LAI = 2");
    _axs[1].plot([0,600], [0,600] .* _max, "k:"; alpha = 0.3);

    _es,_ps = cnpp_light_curve("nature", 4);
    _max = maximum(_ps ./ _es);
    _axs[1].plot(_es, _ps, "k--"; label = "LAI = 4");
    _axs[1].plot([0,600], [0,600] .* _max, "k--"; alpha = 0.3);

    _es,_ps = cnpp_light_curve("nature", 6);
    _max = maximum(_ps ./ _es);
    _axs[1].plot(_es, _ps, "k-"; label = "LAI = 6");
    _axs[1].plot([0,600], [0,600] .* _max, "k-"; alpha = 0.3);
    _axs[1].legend(loc = "upper left");

    # plot optimal light intensity on axis 3
    _es,_ps = cnpp_light_curve("green");
    _max = maximum(_ps ./ _es);
    _axs[4].plot(_es, _ps, "g-");
    _axs[4].plot([0,600], [0,600] .* _max, "g:");

    _es,_ps = cnpp_light_curve("red");
    _max = maximum(_ps ./ _es);
    _axs[4].plot(_es, _ps, "r-");
    _axs[4].plot([0,600], [0,600] .* _max, "r:");

    _es,_ps = cnpp_light_curve("nir");
    _max = maximum(_ps ./ _es);
    _axs[4].plot(_es, _ps, "m-");
    _axs[4].plot([0,600], [0,600] .* _max, "m:");

    # plot optimal light intensity on axis 2
    _lais = collect(1.0:0.5:6.0);
    _maxs = optimal_light_intensity.("nature", _lais);
    _max_lues = [_max[1] for _max in _maxs];
    _max_rads = [_max[2] for _max in _maxs];
    _axs[2].plot(_lais, _max_lues, "k-"; label = "Natural light");
    _axs[3].plot(_lais, _max_rads, "k-"; label = "Natural light");
    _axs[3].legend(loc = "lower right");

    # plot optimal light intensity on axis 4
    _maxs = optimal_light_intensity.("green", _lais);
    _max_lues = [_max[1] for _max in _maxs];
    _max_rads = [_max[2] for _max in _maxs];
    _axs[5].plot(_lais, _max_lues, "g-"; label = "Green light");
    _axs[6].plot(_lais, _max_rads, "g-"; label = "Green light");
    _maxs = optimal_light_intensity.("red", _lais);
    _max_lues = [_max[1] for _max in _maxs];
    _max_rads = [_max[2] for _max in _maxs];
    _axs[5].plot(_lais, _max_lues, "r-"; label = "Red light");
    _axs[6].plot(_lais, _max_rads, "r-"; label = "Red light");
    _maxs = optimal_light_intensity.("nir", _lais);
    _max_lues = [_max[1] for _max in _maxs];
    _max_rads = [_max[2] for _max in _maxs];
    _axs[5].plot(_lais, _max_lues, "m-"; label = "Red+NIR light");
    _axs[6].plot(_lais, _max_rads, "m-"; label = "Red+NIR light");
    _axs[6].legend(loc = "lower right");

    decorate!(_axs;
              add_title = true,
              title_capital = true,
              title_parentheses = false,
              xaxis_labels = ["Light intensity (W m⁻²)", "LAI (-)", "LAI (-)", "Light intensity (W m⁻²)", "LAI (-)", "LAI (-)"],
              xaxis_lims = [(0,600), (1,6), (1,6), (0,600), (1,6), (1,6)],
              yaxis_labels = ["Yield (μmol m⁻² s⁻¹)", "Max efficiency (μmol J⁻¹)", "Opt intensity (W m⁻²)", "Yield (μmol m⁻² s⁻¹)", "Max efficiency (μmol J⁻¹)", "Opt intensity (W m⁻²)"],
              yaxis_lims = [(0,40), (0,0.06), (0,300), (0,40), (0,0.2), (0,150)]);
    save_canvas!(_fig, "6-light"; folder = OUTPUT_FOLDER, formats = ["pdf","png"])

    return nothing
end
