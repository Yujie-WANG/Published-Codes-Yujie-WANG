using Statistics: mean, std

using GriddingMachine.Indexer: lat_ind, lon_ind

using JuliaUtilities.EmeraldRegression: linear_regress
using JuliaUtilities.NetcdfIO: read_nc
using JuliaUtilities.PlotPlants: canvas, decorate!, save_canvas!


CLIMA_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
CLIMA_TAG = "a6_gm1_wd1";
FIG_FOLDER = "";
UNIT_SIF = "mW m⁻² nm⁻¹ sr⁻¹";


function plot_diurnal!()
    _nc_file = "$(CLIMA_FOLDER)/$(CLIMA_TAG)/$(CLIMA_TAG)_2019_1X_1H.c3.rpar.nc";
    _lati = lat_ind(-16.4980; res = 1.0);
    _loni = lon_ind(-56.4120; res = 1.0);
    _inds = read_nc(_nc_file, "ind");
    _fdoy = (_inds .- 0.5 .+ (-56.4120 / 15)) ./ 24;
    _gpps = read_nc(_nc_file, "GPP"   , _loni, _lati);
    _sifs = read_nc(_nc_file, "SIF740", _loni, _lati);
    _pars = read_nc(_nc_file, "PAR"   , _loni, _lati);
    _nirv = read_nc(_nc_file, "NIRvI" , _loni, _lati);
    _nrvp = _nirv .* _pars * 1e6;

    # partition the canvas to 1 big axis and 8 small ones
    _fig = canvas("5_diurnal"; ncol = 4);
    _axs = _fig.axes;

    # plot the panels
    _mask_1 = 100 .<= _fdoy .<= 101;
    _axs[1].plot(0.5:1:24, _gpps[_mask_1]      , "c-"; alpha = 1.0, label = "GPP"  );
    _axs[1].plot(0.5:1:24, _sifs[_mask_1] .* 10, "r-"; alpha = 0.7, label = "SIF"  );
    _axs[1].plot(0.5:1:24, _nirv[_mask_1] .* 30, "b-"; alpha = 0.7, label = "NIRv" );
    _axs[1].plot(0.5:1:24, _nrvp[_mask_1] ./ 25, "b:"; alpha = 0.7, label = "NIRvP");
    _axs[1].legend(loc = "upper left");

    _mask_2 = 321 .<= _fdoy .<= 322;
    _axs[2].plot(0.5:1:24, _gpps[_mask_2]      , "c-"; alpha = 1.0);
    _axs[2].plot(0.5:1:24, _sifs[_mask_2] .* 10, "r-"; alpha = 0.7);
    _axs[2].plot(0.5:1:24, _nirv[_mask_2] .* 30, "b-"; alpha = 0.7);
    _axs[2].plot(0.5:1:24, _nrvp[_mask_2] ./ 25, "b:"; alpha = 0.7);

    _mask_3 = _pars[_mask_1] .> 0;
    _mask_4 = _pars[_mask_2] .> 0;
    _axs[3].plot(_sifs[_mask_1][_mask_3] .* 10, _gpps[_mask_1][_mask_3], "ro"; alpha = 0.7, label = "SIF"  );
    _axs[3].plot(_nrvp[_mask_1][_mask_3] ./ 25, _gpps[_mask_1][_mask_3], "bo"; alpha = 0.7, label = "NIRvP");
    _axs[3].plot(_sifs[_mask_2][_mask_4] .* 10, _gpps[_mask_2][_mask_4], "ro"; alpha = 0.7, mfc = "none");
    _axs[3].plot(_nrvp[_mask_2][_mask_4] ./ 25, _gpps[_mask_2][_mask_4], "bo"; alpha = 0.7, mfc = "none");
    _axs[3].legend(loc = "upper right");

    # create data to plot on the last axis
    _R²_SIF  = Float64[];
    _R²_NIRv = Float64[];
    for _day in 1:365
        _mask_d = _day .<= _fdoy .<= _day + 1;
        _mask_p = _pars[_mask_d] .> 0;
        if sum(_mask_p) >= 4;
            _lr_sif = linear_regress((_sifs[_mask_d][_mask_p],1), _gpps[_mask_d][_mask_p]);
            _lr_nrv = linear_regress((_nrvp[_mask_d][_mask_p],1), _gpps[_mask_d][_mask_p]);
            push!(_R²_SIF , _lr_sif.R²);
            push!(_R²_NIRv, _lr_nrv.R²);
        end;
    end;
    _axs[4].bar([1,2], [mean(_R²_SIF), mean(_R²_NIRv)]; yerr = [std(_R²_SIF), std(_R²_NIRv)], color = ["r", "b"], alpha = 0.7, ecolor = "k", capsize = 10);
    _axs[4].set_xticks([1,2]);
    _axs[4].set_xticklabels(["SIF", "NIRvP"]);

    decorate!(_axs;
              title_labels = String["Wet", "Dry", "", ""],
              xaxis_labels = String["Time of Day 100", "Time of Day 321", "Scaled value", ""],
              yaxis_labels = String["Scaled value", "Scaled value", "GPP (μmol m⁻² s⁻¹)", "R² of GPP ~ X (-)"]);
    decorate!(_axs[1:2]; xaxis_lims = (0,24), xaxis_ticks = collect(0:4:24));
    decorate!(_axs[1:3], yaxis_lims = (-1,38));
    save_canvas!(_fig, "2021_clima_land_fig_5_diurnal"; folder = FIG_FOLDER);

    return nothing
end
