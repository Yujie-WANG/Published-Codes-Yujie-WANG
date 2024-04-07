###############################################################################
#
# Match the soil type and soil water content and leaf water potential
#
###############################################################################
"""
There are additional predawn leaf water potential observations for Ozark flux
    tower site. Thus, we may invert the soil moisture retention curve. This
    function is meant to match the time series of water potential measurements
    with the soil moisture series. Then a fitting is done to generate the best
    fitting retention curve for the soil at Ozark.

    match_ST!(proj::ForestTower2021{FT}) where {FT<:AbstractFloat}

Match the time series of leaf water potential and fit soil van Genuchaten
    parameters, given
- `proj` [`ForestTower2021`](@ref) type project identifier
"""
function match_ST!(proj::ForestTower2021{FT}) where {FT<:AbstractFloat}
    @info tinfo("Reading Flux Tower data and predawn potential data...");
    _ψ_file = "/home/wyujie/RAID/Data/FLUXNET2015/PLWP_MOz.csv";
    _V_file = "/home/wyujie/RAID/Data/FLUXNET2015/US_MOz.csv";
    _ψ_data = read_csv(_ψ_file);
    _V_data = read_csv(_V_file; skiprows=2);
    _ψ_data[!,"SWC"]  = _ψ_data.PLWP .* 0.0;
    _V_data[!,"Year"] = Int.(floor.(_V_data.TIMESTAMP_START ./ 1e8));
    _V_data[!,"DOY"]  = parse_timestamp.(_V_data.TIMESTAMP_START;
                                         in_format="YYYYMMDDhhmm",
                                         out_format="DOY");

    # mask the unrealistic numbers to NaN
    _mask = _V_data.SWC_1_1_1 .< 0;
    _V_data.SWC_1_1_1[_mask] .= NaN;

    @show length(_V_data.Year);
    @show length(_ψ_data.Year);

    @info tinfo("Comparing the data and matching soil water contents...");
    @showprogress for _i in eachindex(_ψ_data.SWC)
        _mask = (_V_data.Year .== _ψ_data.Year[_i]) .*
                (_V_data.DOY .== _ψ_data.DOY[_i]);
        _ψ_data.SWC[_i] = nanmean(_V_data.SWC_1_1_1[_mask]) / 100;
    end

    @info tinfo("Choose !NaN data and fit soil Van Genuchten parameters...");
    _mask  = .!isnan.(_ψ_data.SWC) .* .!isnan.(_ψ_data.PLWP);
    _ψ_obs = _ψ_data.PLWP[_mask];
    _Θ_obs = _ψ_data.SWC[_mask];

    @inline f(x) = (
        _vg = VanGenuchten{Float64}(stype = "Ozark",
                                        α = x[1],
                                        n = x[2],
                                       Θs = 0.45,
                                       Θr = 0.067);
        _ψ_mod = soil_p_25_swc.([_vg], _Θ_obs);
        return -sum( (_ψ_obs .- _ψ_mod) .^ 2 )
    )

    _st = SolutionToleranceND{Float64}([1e-3,1e-4], 30);
    _ms = ReduceStepMethodND{Float64}(x_mins = [1, 1],
                                      x_maxs = [1e4, 10],
                                      x_inis = [600, 1.48],
                                      Δ_inis = [100, 0.1]);
    _αn = find_peak(f, _ms, _st);

    _vg = VanGenuchten{Float64}(stype = "Ozark",
                                    α = _αn[1],
                                    n = _αn[2],
                                   Θs = 0.45,
                                   Θr = 0.067);
    _ψ_mod = soil_p_25_swc.([_vg], _Θ_obs);

    _vg_def = create_soil_VC(VanGenuchten{Float64}(), "Silt Loam");
    _ψ_def  = soil_p_25_swc.([_vg_def], _Θ_obs);

    # plot the results
    _fig,_axs = create_canvas("test", figsize=(10,4));
    _ax1 = _axs[1];
    _ax1.plot(_Θ_obs, _ψ_obs, "+", color="gray");
    _ax1.plot(_Θ_obs, _ψ_mod, ".", color="c");
    _ax1.plot(_Θ_obs, _ψ_def, ".", color="royalblue");

    set_xylabels!(_axs, "Soil water content", "Water potential (MPa)");
    save_canvas!(_fig, "soil_retension_curve.pdf",
                 true);

    return nothing
end








###############################################################################
#
# Match the TROPOMI SIF data with flux tower climate drivers
#
###############################################################################
"""
This function matches the TROPOMI SIF and flux tower data, so that we make 1:1
    comparison.

    match_TROPOMI(
                proj::ForestTower2021{FT},
                site::String,
                year::Int
    ) where {FT<:AbstractFloat}

Align TROPOMI SIF data and flux tower data, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `site` Which site data to match, must be in `NiwotRidge` and `Ozark`
- `year` Which year data to match, best in `2018` and `2019`
"""
function match_TROPOMI(
            proj::ForestTower2021{FT},
            site::String,
            year::Int
) where {FT<:AbstractFloat}
    _nfile = "/home/wyujie/RAID/Data/FLUXNET2015/tropomi/$(site)_$(year).csv";

    # if file exists already, read it directly
    if isfile(_nfile) return read_csv(_nfile) end

    @assert site in ["NiwotRidge", "Ozark"];
    if site == "NiwotRidge"
        _file = "/home/wyujie/RAID/Data/FLUXNET2015/tropomi/US-NR1.csv";
        _lat,_lon = FT(40.03),FT(-105.55);
    else
        _file = "/home/wyujie/RAID/Data/FLUXNET2015/tropomi/US-MOz.csv";
        _lat,_lon = FT(38.74),FT(-92.2);
    end

    # read data
    _data_tower = query_data(proj, year, site, false, "osm");
    _data_tropo = read_csv(_file);

    # read MODIS LAI
    _lai_modis = load_LUT(LAIMODISv006{FT}(), year, "20X", "8D");
    _lais = [read_LUT(_lai_modis,_lat,_lon,i) for i in 1:46];
    _lais = repeat(_lais; inner=8);

    # iterate through the data and save it to a Vector
    _vecs = [];
    _head = [names(_data_tower); "obsSIF"; "Coverage"; "SZA"; "VZA"; "RAA";
             "LAI"];
    for _i in eachindex(_data_tower.TIME)
        # filter the data
        _lsd  = parse_timestamp(_data_tower.TIME[_i]; in_format="YYYYMMDDhhmm",
                                out_format="DATE");
        _lst1 = parse_timestamp(_data_tower.TIME[_i]; in_format="YYYYMMDDhhmm",
                                out_format="FDOY") % 1 * 24;
        _lst2 = _lst1 + 0.5;
        _doyn = parse_timestamp(_data_tower.TIME[_i]; in_format="YYYYMMDDhhmm",
                                out_format="DOY");
        _sel1 = (_data_tropo.DATE .== _lsd) .*
                (_lst1 .<= _data_tropo.LST .<= _lst2) .*
                (_data_tropo.CF .< 0.1);

        # save the data
        for _j in 1:sum(_sel1)
            _data_sel = _data_tropo[_sel1,:];
            _vec = [_data_tower[_i,:]...; _data_sel.SIF[_j]; _data_sel.PER[_j];
                    _data_sel.SZA[_j]; _data_sel.VZA[_j]; _data_sel.RAA[_j];
                    _lais[_doyn]];
            push!(_vecs, _vec);
        end
    end

    # save the new file as a CSV
    save_csv!(_nfile, _vecs, _head);

    return read_csv(_nfile)
end
