using PyCall: pyimport

using GriddingMachine.Indexer: lat_ind, lon_ind

using JuliaUtilities.NetcdfIO: read_nc
using JuliaUtilities.PlotPlants: canvas, decorate!, save_canvas!

CCRS = pyimport("cartopy.crs");


CLIMA_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
CLIMA_TAG = "a6_gm1_wd1";
FIG_FOLDER = "";
UNIT_SIF = "mW m⁻² nm⁻¹ sr⁻¹";


function plot_examples!()
    _nc_file = "$(CLIMA_FOLDER)/$(CLIMA_TAG)/$(CLIMA_TAG)_2019_1X_1D.c3.rpar.nc";

    # choose 8 flux tower sites
    _lats = [  40.0329, -16.4980,  5.2685, 51.0792, -25.0197,  62.2550,  26.7414, -35.6566];
    _lons = [-105.5464, -56.4120, -2.6942, 10.4522,  31.4969, 129.1680, 115.0581, 148.1517];
    _labs = ["a", "b", "c", "d", "e", "f", "g", "h"];
    _tits = ["US-NR1", "BR-Npw", "GH-Ank", "DE-Hai", "ZA-Kru", "RU-Skp", "CN-Qia", "AU-Tum"];
    _lati = lat_ind.(_lats; res = 1.0);
    _loni = lon_ind.(_lons; res = 1.0);

    # partition the canvas to 1 big axis and 8 small ones
    _fig = canvas("2_example"; nrow = 4, ncol = 4, axids = [3,4,7,8,9,10,11,12], figsize = (15,9));
    _axs = _fig.axes;
    _txs = [_ax.twinx() for _ax in _axs];
    _map = _fig.add_subplot(2,2,1; projection = CCRS.PlateCarree());

    for _i in eachindex(_lats)
        _map.plot(_lons[_i], _lats[_i], "ro", mfc = "none", markersize = 12, transform = CCRS.PlateCarree());
        _map.text(_lons[_i], _lats[_i], _labs[_i], ha = "center", va = "center", color = "r", transform = CCRS.PlateCarree());
        _gpp = read_nc(_nc_file, "mGPP"   , _loni[_i], _lati[_i]);
        _sif = read_nc(_nc_file, "mSIF740", _loni[_i], _lati[_i]);
        _niv = read_nc(_nc_file, "mNIRvI" , _loni[_i], _lati[_i]);
        _axs[_i].plot(1:365, _gpp    , "c-", alpha = 1.0, label = "GPP");
        _txs[_i].plot(1:365, _sif    , "r-", alpha = 0.5, label = "SIF at 740 nm");
        _txs[_i].plot(1:365, _niv * 2, "b-", alpha = 0.5, label = "2×NIRv");
    end;
    _axs[1].legend(loc = "upper left");
    _txs[1].legend(loc = "upper right");
    _map.stock_img();
    _map.set_global();

    decorate!(_axs; title_labels = _tits, xaxis_lims = (0,367), yaxis_lims = (0,20));
    decorate!(_txs; yaxis_lims = (0,2));
    decorate!(_axs[5:8]; xaxis_labels = "Day of year 2019");
    decorate!(_axs[[1,3,5]]; yaxis_labels = "GPP (g C m⁻² d⁻¹)");
    decorate!(_txs[4:4]; yaxis_labels = "NIRv (-) or SIF ($(UNIT_SIF))");
    save_canvas!(_fig, "2021_clima_land_fig_2_example"; folder = FIG_FOLDER);

    return nothing
end;
