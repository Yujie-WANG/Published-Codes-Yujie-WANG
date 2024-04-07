using CairoMakie: save
using GeoMakie: Figure, GeoAxis

using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: read_LUT

using JuliaUtilities.NetcdfIO: read_nc


include("plot_axis.jl");


CLIMA_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
CLIMA_TAG = "a6_gm1_wd1";
FIG_FOLDER = "";
UNIT_A = "μmol CO₂ m⁻² s⁻¹";
UNIT_SIF = "mW m⁻² sr⁻¹ nm⁻¹";


function plot_maps!()
    _1m_file = "$(CLIMA_FOLDER)/$(CLIMA_TAG)/$(CLIMA_TAG)_2019_1X_1M.c3.rpar.nc";
    _8d_file = "$(CLIMA_FOLDER)/$(CLIMA_TAG)/$(CLIMA_TAG)_2019_1X_8D.c3.rpar.nc";

    _fig = Figure(resolution=(1800,1200), fontsize=24);
    _ax11 = GeoAxis(_fig[1,1]);
    _cx11 = _fig[1,2];
    _ax12 = GeoAxis(_fig[1,3]);
    _cx12 = _fig[1,4];
    _ax21 = GeoAxis(_fig[2,1]);
    _cx21 = _fig[2,2];
    _ax22 = GeoAxis(_fig[2,3]);
    _cx22 = _fig[2,4];
    _ax31 = GeoAxis(_fig[3,1]);
    _cx31 = _fig[3,2];
    _ax32 = GeoAxis(_fig[3,3]);
    _cx32 = _fig[3,4];

    _lats = read_nc(_1m_file, "lat");
    _lons = read_nc(_1m_file, "lon");

    # 1. plot the GPP comparison
    _ref_2x,_ = read_LUT(query_collection("GPP_MPI_RS_2X_1M_2019_V1"));
    _ref_2x ./= (1e-6 * 3600 * 24 * 12);
    _ref_2x[isnan.(_ref_2x)] .= 0;
    _ref_1x = regrid(_ref_2x, 1);
    _clm_1x = read_nc(_1m_file, "mGPP");
    _clm_1x[isnan.(_clm_1x)] .= 0;
    _ref_map,_clm_map = processed_data_2d(_ref_1x, _clm_1x);
    plot_map_on_axis!(_ax11, _cx11, _lons, _lats, _ref_map, "MPI RS GPP\n($(UNIT_A))", vminvmax=(0,8));
    plot_map_on_axis!(_ax12, _cx12, _lons, _lats, _clm_map, "CliMA GPP\n($(UNIT_A))", vminvmax=(0,15));

    # 2. plot the SIF 740 comparison
    _ref_1x,_ = read_LUT(query_collection("SIF_TROPOMI_740DC_1X_1M_2019_V1"));
    _ref_1x[isnan.(_ref_1x)] .= 0;
    _clm_1x = read_nc(_1m_file, "mSIF740");
    _clm_1x[isnan.(_clm_1x)] .= 0;
    _ref_map,_clm_map = processed_data_2d(_ref_1x, _clm_1x);
    plot_map_on_axis!(_ax21, _cx21, _lons, _lats, _ref_map, "TROPOMI SIF₇₄₀\n($(UNIT_SIF))", vminvmax=(0,0.63));
    plot_map_on_axis!(_ax22, _cx22, _lons, _lats, _clm_map, "CliMA SIF₇₄₀\n($(UNIT_SIF))", vminvmax=(0,0.85));

    # 3. plot the NIRv comparison
    _ref_1x = read_nc("/net/fluo/data1/data/ERA5/simulation/modis_mcd43_1x_8d_2019.nc", "NIRv");
    _ref_1x[isnan.(_ref_1x)] .= 0;
    _ref_1x .*= LM_1X;
    _clm_1x = read_nc(_8d_file, "mNIRvI");
    _clm_1x[isnan.(_clm_1x)] .= 0;
    _ref_map,_clm_map = processed_data_2d(_ref_1x, _clm_1x);
    plot_map_on_axis!(_ax31, _cx31, _lons, _lats, _ref_map, "MODIS NIRv (-)", vminvmax=(0,0.4));
    plot_map_on_axis!(_ax32, _cx32, _lons, _lats, _clm_map, "CliMA NIRv (-)", vminvmax=(0,0.605));

    # save("$(FIG_FOLDER)/3_maps.pdf", _fig)
    save("$(FIG_FOLDER)/2021_clima_land_fig_3_maps.png", _fig);

    return nothing
end;
