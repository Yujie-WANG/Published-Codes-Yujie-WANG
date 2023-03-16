using JuliaUtilities.NetcdfIO: read_nc
using JuliaUtilities.PlotPlants: canvas, save_canvas!

using Land.EmeraldConstants: M_H₂O
using Land.WaterPhysics: latent_heat_vapor


include("plot_axis.jl");


CLIMA_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
CLIMA_TAG = "a6_gm1_wd1";
FIG_FOLDER = "";
UNIT_A = "μmol CO₂ m⁻² s⁻¹";
UNIT_LE = "W m⁻²";
UNIT_SIF = "mW m⁻² sr⁻¹ nm⁻¹";


function plot_hexbins!()
    _1m_2010 = "$(CLIMA_FOLDER)/$(CLIMA_TAG)/$(CLIMA_TAG)_2010_1X_1M.c3.rpar.nc";
    _1m_2019 = "$(CLIMA_FOLDER)/$(CLIMA_TAG)/$(CLIMA_TAG)_2019_1X_1M.c3.rpar.nc";
    _8d_2019 = "$(CLIMA_FOLDER)/$(CLIMA_TAG)/$(CLIMA_TAG)_2019_1X_8D.c3.rpar.nc";
    _ms_file = "/net/fluo/data1/data/ERA5/simulation/modis_mcd43_1x_8d_2019.nc";

    _fig = canvas("hexbin", nrow=3, ncol=3, figsize=(12,9.5));
    _axs = _fig.axes;

    # 1. plot the GPP comparison
    _ref_2x,_ = read_LUT(query_collection("GPP_MPI_RS_2X_1M_2019_V1"));
    _ref_2x ./= (1e-6 * 3600 * 24 * 12);
    _ref_2x[isnan.(_ref_2x)] .= 0;
    _ref_1x = regrid(_ref_2x, 1);
    _clm_1x = read_nc(_1m_2019, "mGPP");
    _clm_1x[isnan.(_clm_1x)] .= 0;

    plot_comparison!(_fig, _axs[1], _ref_1x, _clm_1x);
    _axs[1].set_xlabel("MPI RS GPP ($(UNIT_A))");
    _axs[1].set_ylabel("CliMA GPP ($(UNIT_A))");

    # 2. plot the LE comparison
    _ref_2x,_ = read_LUT(query_collection("LE_MPI_RS_2X_1M_2010_V1"));
    _ref_2x[isnan.(_ref_2x)] .= 0;
    _ref_1x = regrid(_ref_2x, 1);
    _clm_1x = read_nc(_1m_2010, "mT");
    _clm_1x[isnan.(_clm_1x)] .= 0;
    _ref_1x .*= 1e6 / 24 / 3600;
    _clm_1x .*= latent_heat_vapor(298.15) * M_H₂O();

    plot_comparison!(_fig, _axs[2], _ref_1x, _clm_1x);
    _axs[2].set_xlabel("MPI RS LE ($(UNIT_LE) )");
    _axs[2].set_ylabel("CliMA T LE ($(UNIT_LE))");

    # 3. plot the SIF 683 comparison
    _ref_5x,_ = read_LUT(query_collection("SIF_TROPOMI_683DC_5X_1M_2019_V2"));
    _ref_5x[isnan.(_ref_5x)] .= 0;
    _ref_1x = regrid(_ref_5x, 1);
    _clm_1x = read_nc(_1m_2019, "mSIF683");
    _clm_1x[isnan.(_clm_1x)] .= 0;

    plot_comparison!(_fig, _axs[3], _ref_1x, _clm_1x);
    _axs[3].set_xlabel("TROPOMI SIF\$_{683}\$ ($(UNIT_SIF))");
    _axs[3].set_ylabel("CliMA SIF\$_{683}\$ ($(UNIT_SIF))");

    # 4. plot the SIF 740 comparison
    _ref_1x,_ = read_LUT(query_collection("SIF_TROPOMI_740DC_1X_1M_2019_V1"));
    _ref_1x[isnan.(_ref_1x)] .= 0;
    _clm_1x = read_nc(_1m_2019, "mSIF740");
    _clm_1x[isnan.(_clm_1x)] .= 0;

    plot_comparison!(_fig, _axs[4], _ref_1x, _clm_1x);
    _axs[4].set_xlabel("TROPOMI SIF\$_{740}\$ ($(UNIT_SIF))");
    _axs[4].set_ylabel("CliMA SIF\$_{740}\$ ($(UNIT_SIF))");

    # 5. plot the SIF 757 comparison
    _ref_5x,_ = read_LUT(query_collection("SIF_OCO2_757DC_5X_1M_2019_V3"));
    _ref_1x = regrid(_ref_5x, 1);
    _ref_1x[isnan.(_ref_1x)] .= 0;
    _clm_1x = read_nc(_1m_2019, "mSIF757");
    _clm_1x[isnan.(_clm_1x)] .= 0;

    plot_comparison!(_fig, _axs[5], _ref_1x, _clm_1x);
    _axs[5].set_xlabel("OCO-2 SIF\$_{757}\$ ($(UNIT_SIF))");
    _axs[5].set_ylabel("CliMA SIF\$_{757}\$ ($(UNIT_SIF))");

    # 6. plot the SIF 771 comparison
    _ref_5x,_ = read_LUT(query_collection("SIF_OCO2_771DC_5X_1M_2019_V3"));
    _ref_1x = regrid(_ref_5x, 1);
    _ref_1x[isnan.(_ref_1x)] .= 0;
    _clm_1x = read_nc(_1m_2019, "mSIF771");
    _clm_1x[isnan.(_clm_1x)] .= 0;

    plot_comparison!(_fig, _axs[6], _ref_1x, _clm_1x);
    _axs[6].set_xlabel("OCO-2 SIF\$_{771}\$ ($(UNIT_SIF))");
    _axs[6].set_ylabel("CliMA SIF\$_{771}\$ ($(UNIT_SIF))");

    # 7. plot the NDVI comparison
    _ref_1x = read_nc(_ms_file, "NDVI");
    _ref_1x[isnan.(_ref_1x)] .= 0;
    _clm_1x = read_nc(_8d_2019, "mNDVI");
    _clm_1x[isnan.(_clm_1x)] .= 0;

    plot_comparison!(_fig, _axs[7], _ref_1x, _clm_1x);
    _axs[7].set_xlabel("MODIS NDVI (-)");
    _axs[7].set_ylabel("CliMA NDVI (-)");

    # 8. plot the EVI comparison
    _ref_1x = read_nc(_ms_file, "EVI");
    _ref_1x[isnan.(_ref_1x)] .= 0;
    _clm_1x = read_nc(_8d_2019, "mEVI");
    _clm_1x[isnan.(_clm_1x)] .= 0;

    plot_comparison!(_fig, _axs[8], _ref_1x, _clm_1x);
    _axs[8].set_xlabel("MODIS EVI (-)");
    _axs[8].set_ylabel("CliMA EVI (-)");

    # 9. plot the NIRv comparison
    _ref_1x = read_nc(_ms_file, "NIRv");
    _ref_1x[isnan.(_ref_1x)] .= 0;
    _ref_1x .*= LM_1X;
    _clm_1x = read_nc(_8d_2019, "mNIRvI");
    _clm_1x[isnan.(_clm_1x)] .= 0;

    plot_comparison!(_fig, _axs[9], _ref_1x, _clm_1x);
    _axs[9].set_xlabel("MODIS NIRv (-)");
    _axs[9].set_ylabel("CliMA NIRv (-)");

    # save the figure
    save_canvas!(_fig, "2021_clima_land_fig_4_hexbins"; folder = FIG_FOLDER);

    return nothing
end;
