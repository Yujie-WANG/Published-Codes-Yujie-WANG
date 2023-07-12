using Emerald.EmeraldLand.LeafOptics: leaf_spectra!
using Emerald.EmeraldLand.Namespace: HyperspectralAbsorption, HyperspectralLeafBiophysics, MultiLayerSPAC, SPACConfiguration, WaveLengthSet
using Emerald.EmeraldLand.SPAC: GPP, PPAR, initialize!, soil_plant_air_continuum!, update!
using Emerald.EmeraldIO.Text: read_csv

using Emerald.EmeraldVisualization: canvas, decorate!, save_canvas!


# mkpath if the path does not exist
if !isdir("$(@__DIR__)/output")
    mkpath("$(@__DIR__)/output");
end;
OUTPUT_FOLDER = "$(@__DIR__)/output";


# function to plot spectra for A, B, D, and F
function plot_spectra!()
    # read the data from CSV file
    _data = read_csv("$(@__DIR__)/../../data/spectra.csv");

    # create a canvas to plot the spectra
    _fig = canvas("2"; dpi = 600, figsize = (5,3));
    _axs = _fig.axes;

    _chls = ["A", "B", "D", "F"];
    _cols = ["orange", "skyblue", "salmon", "red"];
    for _i in eachindex(_chls)
        _xs = _data[:,"WL$(_chls[_i])"];
        _ys = _data[:,"CHL$(_chls[_i])"];
        _xx = _xs[.!ismissing.(_xs)];
        _yy = _ys[.!ismissing.(_ys)];
        _axs[1].plot(_xx, _yy, "-"; color = _cols[_i], label = lowercase("chl $(_chls[_i])"));
    end;
    _axs[1].legend(loc = "upper left", ncol = 2);

    decorate!(_axs; add_title = false, xaxis_lims = (370,800), yaxis_lims = (0,1.1), yaxis_ticks = []);
    save_canvas!(_fig, "2-spectra"; folder = OUTPUT_FOLDER, formats = ["pdf","png"]);

    return nothing
end


# function to generate artificial absorption curve
function artificial_curve()
    # generate a truth curve
    _config = SPACConfiguration{Float64}();
    _wls = _config.WLSET;
    _bha = HyperspectralAbsorption{Float64}();
    _bio = HyperspectralLeafBiophysics(_config);
    _mha = HyperspectralAbsorption{Float64}();
    _mio = HyperspectralLeafBiophysics(_config);

    # change the shape of the absorption feature of chl
    for _i in eachindex(_wls.Λ)[end:-1:1]
        if 680 < _wls.Λ[_i] < 775
            _mha.K_CAB[_i] = _mha.K_CAB[_i-5];
        end;
    end;

    leaf_spectra!(_bio, _wls, _bha, Float64(50));
    leaf_spectra!(_mio, _wls, _mha, Float64(50));

    return _wls.Λ, _bio.α_sw, _mio.α_sw
end


# function to plot spectra for wild type and mutant
function plot_absorption!()
    # read the data from CSV file
    (_wls,_bio,_mio) = artificial_curve();

    # create a canvas to plot the spectra
    _fig = canvas("2"; dpi = 600, figsize = (5,3));
    _axs = _fig.axes;

    _axs[1].plot(_wls, _bio; color = "skyblue", label = "wild");
    _axs[1].plot(_wls, _mio; color = "salmon", label = "mutant");
    _axs[1].legend(loc = "lower left", ncol = 2);

    decorate!(_axs; add_title = false, xaxis_lims = (400,800), yaxis_lims = (0,1), yaxis_ticks = []);
    save_canvas!(_fig, "2-absorption"; folder = OUTPUT_FOLDER, formats = ["pdf","png"]);

    return nothing
end


# plot the scenario of uniform absorption coefficient
function plot_uniform_absorption!()
    @inline gpp_func(α) = (
        _config = SPACConfiguration{Float64}();
        _spac = MultiLayerSPAC(_config);
        initialize!(_spac, _config);
        for _leaf in _spac.LEAVES
            leaf_spectra!(_leaf.BIO, _config.WLSET, (1-α)/2, (1-α)/2, (1-α)/2, (1-α)/2);
        end;
        soil_plant_air_continuum!(_spac, _config, 3600.0; p_on = false, t_on = false, θ_on = false);
        soil_plant_air_continuum!(_spac, _config, 1.0; p_on = false, t_on = false, θ_on = false);

        return GPP(_spac)
    );

    _αs = collect(0.02:0.02:0.98);
    _ys = gpp_func.(_αs);
    _mask = (_ys .> maximum(_ys) * 0.99);

    # create a canvas to plot the spectra
    _fig = canvas("2"; dpi = 600, figsize = (7,3));
    _axs = _fig.axes;

    _axs[1].plot(_αs, _ys; color = "black");
    _axs[1].fill_between(_αs[_mask], 0, 32; color = "forestgreen", alpha = 0.3);

    decorate!(_axs; add_title = false, xaxis_labels = "Uniform Absorption Coefficient", xaxis_lims = (0,1), yaxis_labels = "ΣA", yaxis_lims = (0,32));
    save_canvas!(_fig, "2-coefficient"; folder = OUTPUT_FOLDER, formats = ["pdf","png"]);

    return nothing
end


# function to plot spectra for wild type and mutant
function plot_absorption_shaded!()
    # read the data from CSV file
    (_wls,_bio,_mio) = artificial_curve();

    # create a canvas to plot the spectra
    _fig = canvas("2"; dpi = 600, figsize = (7,3));
    _axs = _fig.axes;

    _axs[1].plot(_wls, _bio; color = "skyblue", label = "wild");
    _axs[1].plot(_wls, _mio; color = "salmon", label = "mutant");
    _axs[1].fill_between(_wls, 0.44, 0.64; color = "forestgreen", alpha = 0.3);
    _axs[1].legend(loc = "lower left", ncol = 2);

    decorate!(_axs; add_title = false, xaxis_lims = (400,800), yaxis_lims = (0,1), yaxis_ticks = []);
    save_canvas!(_fig, "2-absorption-shaded"; folder = OUTPUT_FOLDER, formats = ["pdf","png"]);

    return nothing
end


# function to test the multi-farming
function test_nir_mutant!()
    # default setup
    _config = SPACConfiguration{Float64}();
    _spac = MultiLayerSPAC(_config);
    initialize!(_spac, _config);
    update!(_spac, _config; lai = 4, vcmax_expo = 0.15);
    for _ in 1:60
        soil_plant_air_continuum!(_spac, _config, 60.0; p_on = false, t_on = false, θ_on = false);
    end;
    @info "Default" GPP(_spac) PPAR(_spac);

    # lower canopy with NIR mutant
    _config = SPACConfiguration{Float64}();
    _spac = MultiLayerSPAC(_config);
    initialize!(_spac, _config);
    update!(_spac, _config; lai = 4, vcmax_expo = 0.15);
    # change the shape of the absorption feature of chl
    _wls = _config.WLSET;
    _lha = deepcopy(_config.LHA);
    for _i in eachindex(_wls.Λ)[end:-1:1]
        if 680 < _wls.Λ[_i] < 775
            _lha.K_CAB[_i] = _lha.K_CAB[_i-5];
        end;
    end;
    for _i in eachindex(_spac.LEAVES)
        _spac.LEAVES[_i].BIO._v_storage = 0;
        leaf_spectra!(_spac.LEAVES[_i].BIO, _wls, _lha, _spac.LEAVES[_i].HS.v_storage);
    end;
    for _ in 1:60
        soil_plant_air_continuum!(_spac, _config, 60.0; p_on = false, t_on = false, θ_on = false);
    end;
    @info "NIR mutant" GPP(_spac) PPAR(_spac);

    return nothing
end
#=
Info: Default
  GPP(_spac) = 28.279960730392574
  PPAR(_spac) = 1166.7314760019908
Info: NIR mutant
  GPP(_spac) = 28.351008114398894
  PPAR(_spac) = 1230.3745728553647
=#
