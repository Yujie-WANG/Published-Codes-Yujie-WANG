###############################################################################
#
# Query data frame from flux tower dataset
#
###############################################################################
"""
Query data to use with simulation functions. Note that the CSV files are from
    deployed artifact. Remember to update the artifact if you have changed the
    input CSV files.
"""
function query_data end




"""
This method query the data frame required to run the simulations:

    query_data(proj::ForestTower2021{FT},
               year::Int,
               site::String = "NiwotRidge",
               filtering::Bool = true,
               label::String = "osm"
    ) where {FT<:AbstractFloat}

Query DataFrame formated data, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `year` Which year of data to query
- `site` Site of flux tower, must be `NiwotRidge` or `Ozark`
- `filtering` True by default. If true, filtering the growing season only,
- `label` Label to identify the simulation
"""
query_data(proj::ForestTower2021{FT},
           year::Int,
           site::String = "NiwotRidge",
           filtering::Bool = true,
           label::String = "osm"
) where {FT<:AbstractFloat} =
(
    # make sure that site is supported
    @assert site in ["NiwotRidge", "Ozark"];

    # read data from artifact
    _path = artifact"2021_forest_tower_data" * "/$(site)_$(year).csv";
    if !isfile(_path)
        return nothing
    end;

    # continue if file exists
    _df = read_csv(_path);

    # create new columns in the new dataframe
    _df.Day   = parse_timestamp.(_df.TIME; in_format="YYYYMMDDhhmm",
                                 out_format="DOY");
    _df.Hour  = [parse(Int, string(_df.TIME[i])[9:10])
                 for i in eachindex(_df.TIME)];
    _df.Minu  = [parse(Int, string(_df.TIME[i])[11:12]) + 15
                 for i in eachindex(_df.TIME)];
    _df.ObsC  = _df.F_CO2 .* -1;
    _df.ObsE  = _df.LE_H2O ./ latent_heat_vapor(298.15) ./ M_H₂O();
    _df.ModC  = [NaN for i in eachindex(_df.TIME)];
    _df.ModE  = [NaN for i in eachindex(_df.TIME)];
    _df.Site  = [site for i in eachindex(_df.TIME)];
    _df.Label = [label for i in eachindex(_df.TIME)];

    # filter the days only for those with seven day mean >= 1
    if filtering
        _day_start = 1;
        _day_end = 366;
        for _day in 1:366
            _mask = (_day .<= _df.Day .<= _day+7);
            if nanmean(_df.ObsC[_mask]) >= 1
                _day_start = _day;
                break
            end;
        end;
        for _day in 366:-1:1
            _mask = (_day-7 .<= _df.Day .<= _day);
            if nanmean(_df.ObsC[_mask]) >= 1
                _day_end = _day;
                break
            end;
        end;
        _mask = _day_start .< _df.Day .< _day_end;

        return _df[_mask,:]
    else
        return _df
    end;
)




"""
This method query the parameter vector required to run the simulations in
    parallel:

    query_data(proj::ForestTower2021{FT},
               debugging::Bool = true;
               years::Array{Int,1} = collect(2006:2020),
               sites::Vector{String} = ["NiwotRidge", "Ozark"],
               labels::Vector{String} = ["bbm", "bbmg", "med", "medg", "osm"]
    ) where {FT<:AbstractFloat}

Query parameter vector to run the simulations in parallel, given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `debugging` If true, show workers output; if false, show progress bar
- `years` Which years of data to query
- `sites` Sites of flux tower, must be `NiwotRidge` or `Ozark`
- `labels` Labels to identify the simulation

Each item in the parameter vector includes the following
- Project indentifier
- `SPACMono` type struct
- DataFrame of weather data
- Year
- Debugging mode
- Site
- Model label
"""
query_data(proj::ForestTower2021{FT},
           debugging::Bool = true;
           years::Array{Int,1} = collect(2006:2020),
           sites::Vector{String} = ["NiwotRidge", "Ozark"],
           labels::Vector{String} = ["bbm", "bbmg", "med", "medg", "osm"]
) where {FT<:AbstractFloat} =
(
    @info tinfo("Creating parameter vector to pass to threadings...");

    _params = [];
    _osm = OSMWang{FT}();
    _bbm = ESMBallBerry{FT}(g0=0.001, g1=9.0);
    for _site in sites
        if _site == "NiwotRidge"
            _med = ESMMedlyn{FT}(g0=0.001, g1=74.31);
            _chl = FT(46.82);
        elseif _site=="Ozark"
            _med = ESMMedlyn{FT}(g0=0.001, g1=140.72);
            _chl = FT(57.23);
        end;

        _sms = [_bbm, _bbm, _med, _med, _osm];
        for _i in eachindex(labels)
            _label = labels[_i];
            _sm = _sms[_i];
            _node = create_spac(proj, _site, _sm, FT(20), FT(0.3), _chl);
            for _year in years
                _df = query_data(proj, _year, _site, true, _label);
                if _df !== nothing
                    _p = [proj, _node, _df, _year, debugging, _site, _label];
                    push!(_params, deepcopy(_p));
                end;
            end;
        end;
    end;

    return _params
)
