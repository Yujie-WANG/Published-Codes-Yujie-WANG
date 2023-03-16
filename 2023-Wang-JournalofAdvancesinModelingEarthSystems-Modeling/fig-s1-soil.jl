#=
function plot_swc_comparison()
    # lats and lons from the example figure
    _era5_data = read_csv("/home/wyujie/Projects/2021_clima_land/era5_2019_55_329_1X.csv");

    # process flux tower data
    # _flux_data = read_csv("/home/wyujie/RAID/AmeriFlux/FLX_AU-Tum_FLUXNET2015_FULLSET_HR_2001-2014_2-4.csv");
    # _flux_subd = _flux_data[:,["TIMESTAMP_START", "TIMESTAMP_END", "SWC_F_MDS_1", "SWC_F_MDS_1_QC"]];
    # save_csv!(_flux_subd, "/home/wyujie/RAID/Data/LandGPP/AU-Tum.csv");
    # @show names(_era5_data);
    # @show names(_flux_data);
    _flux_data = read_csv("/home/wyujie/Projects/2021_clima_land/AU-Tum.csv");
    _flux_data.DOY  = [parse_timestamp(_flux_data.TIMESTAMP_START[_i]; in_format="YYYYMMDDhhmm", out_format="DOY") for _i in eachindex(_flux_data.TIMESTAMP_START)];
    _flux_data.FDOY = [parse_timestamp(_flux_data.TIMESTAMP_START[_i]; in_format="YYYYMMDDhhmm", out_format="FDOY") for _i in eachindex(_flux_data.TIMESTAMP_START)];
    _flux_data.SWC  = _flux_data.SWC_F_MDS_1 ./ 100;
    _flux_data.SWC[_flux_data.SWC .< 0] .= NaN;
    _flux_data.SWC[_flux_data.SWC .< 0] .= NaN;

    # calculate the daily mean
    _mean_flux = ones(365) .* NaN;
    _stdv_flux = ones(365) .* NaN;
    _mean_era1 = ones(365) .* NaN;
    _mean_era2 = ones(365) .* NaN;
    _mean_era3 = ones(365) .* NaN;
    _mean_era4 = ones(365) .* NaN;
    _mean_era5 = ones(365) .* NaN;
    _stdv_era5 = ones(365) .* NaN;
    for _doy in 1:365
        _mean_flux[_doy] = nanmean(_flux_data.SWC[_flux_data.DOY .== _doy]);
        _stdv_flux[_doy] = nanstd(_flux_data.SWC[_flux_data.DOY .== _doy]);
        _mean_era1[_doy] = nanmean(_era5_data.swvl1[_doy - 1 .<= _era5_data.FDOY .<= _doy]);
        _mean_era2[_doy] = nanmean(_era5_data.swvl2[_doy - 1 .<= _era5_data.FDOY .<= _doy]);
        _mean_era3[_doy] = nanmean(_era5_data.swvl3[_doy - 1 .<= _era5_data.FDOY .<= _doy]);
        _mean_era4[_doy] = nanmean(_era5_data.swvl4[_doy - 1 .<= _era5_data.FDOY .<= _doy]);
        _mean_era5[_doy] = nanmean([_mean_era1[_doy], _mean_era2[_doy], _mean_era3[_doy], _mean_era4[_doy]]);
        _stdv_era5[_doy] = nanstd([_mean_era1[_doy], _mean_era2[_doy], _mean_era3[_doy], _mean_era4[_doy]]);
    end;

    # calculate the soil water potential based on the soil type (commented command below are based on Griddingmachine v0.1)
    # _s_α  = load_LUT(VGMAlphaJules{Float64}(), "12X", "1Y", 1);
    # _s_ln = load_LUT(VGMLogNJules{Float64}(), "12X", "1Y", 1);
    # _s_Θr = load_LUT(VGMThetaRJules{Float64}(), "12X", "1Y", 1);
    # _s_Θs = load_LUT(VGMThetaSJules{Float64}(), "12X", "1Y", 1);
    # _α  = read_LUT(_s_α, -35.6566, 148.1517, 1);
    # _n  = exp( read_LUT(_s_ln, -35.6566, 148.1517, 1) )
    # _Θr = read_LUT(_s_Θr, -35.6566, 148.1517, 1);
    # _Θs = read_LUT(_s_Θs, -35.6566, 148.1517, 1);
    _sh = VanGenuchten{Float64}(stype = "AU-Tum", α = 480.19061226314966, n = 1.1282569798047912, Θs = 0.404426842307051, Θr = 0.07456589782507056);

    # plot the time series of SWC and Psoil
    _fig,_axs = create_canvas(FIG_SWC; nrow=2, figsize=(7,6.5));
    _ax1,_ax2 = _axs;

    _ax1.plot(1:365, _mean_era5, "k-"; label="ERA5", linewidth=2);
    #_ax1.fill_between(1:365, _mean_era5 .+ _stdv_era5, _mean_era5 .- _stdv_era5, color="k", alpha=0.2);
    _ax1.plot(1:365, _mean_flux, "r-"; label="Flux tower", linewidth=2);
    _ax1.fill_between(1:365, _mean_flux .+ _stdv_flux, _mean_flux .- _stdv_flux, color="r", alpha=0.2);
    _ax1.legend(loc="upper right", ncol=2);

    _ax2.plot(1:365, soil_p_25_swc.([_sh], _mean_era5), "k-"; linewidth=2);
    #_ax2.fill_between(1:365, soil_p_25_swc.([_sh], _mean_era5 .+ _stdv_era5), soil_p_25_swc.([_sh], _mean_era5 .- _stdv_era5), color="k", alpha=0.2);
    _ax2.plot(1:365, soil_p_25_swc.([_sh], _mean_flux), "r-"; linewidth=2);
    _ax2.fill_between(1:365, soil_p_25_swc.([_sh], _mean_flux .+ _stdv_flux), soil_p_25_swc.([_sh], _mean_flux .- _stdv_flux), color="r", alpha=0.2);

    set_xylabels!(_axs, ["", "Day of year (-)"], ["SWC (-)", "\$\\Psi_\\text{soil}\$ (MPa)"]);
    set_xylims!(_axs, [-1,367], [[0.12, 0.42], [-5, 0]]);
    save_canvas!(_fig, "$(FIG_FOLDER)/$(FIG_SWC).pdf", SAVING);

    return _fig
end;
=#
