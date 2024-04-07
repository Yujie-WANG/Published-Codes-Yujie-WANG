###############################################################################
#
# Test the impact of atmospheric pressure on carbon and water flux
#
###############################################################################
"""
This function tests how much atmospheric pressure may impact the carbon and
    water flux. This function compares two scenarios for Niwot Ridge flux tower
    site, and the simulation results are then saved to `df_sea_level.csv` and
    `df_flux_tower.csv`:
- Mean sea level atmospheric pressure is used
- Mean stand atmospheric pressure is used

    test_patm!(proj::ForestTower2021{FT}) where {FT<:AbstractFloat}

Run annual simulation assuming different `p_atm` scenarios, given
- `proj` [`ForestTower2021`](@ref) type project identifier
"""
function test_patm!(proj::ForestTower2021{FT}) where {FT<:AbstractFloat}
    _df_atm     = query_data(proj, 2019, "NiwotRidge", false, "osm");
    _df_atm.DOY = _df_atm.Day .+ (_df_atm.Hour .+ _df_atm.Minu ./ 60) ./ 24;
    _df_sea     = deepcopy(_df_atm);
    _rbase      = Q10TD{FT}(5.6, 298.15, 1.7);
    _node       = create_spac(proj, "NiwotRidge", OSMWang{FT}(), FT(15),
                              FT(0.13), FT(46.82));
    initialize_spac_canopy!(_node);

    # run simulations
    simulation!(proj, deepcopy(_node), _df_atm, "NiwotRidge", true , _rbase);
    simulation!(proj, deepcopy(_node), _df_sea, "NiwotRidge", false, _rbase);

    # save the simulation to local CSV files
    _folder = "/home/wyujie/RAID/Data/FLUXNET2015/simulation";
    save_csv!("$(_folder)/df_flux_tower.csv", _df_atm);
    save_csv!("$(_folder)/df_sea_level.csv" , _df_sea);

    return nothing;
end
