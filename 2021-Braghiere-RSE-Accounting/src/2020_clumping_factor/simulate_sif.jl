###############################################################################
#
# Run the SIF simulation for Renato
#
###############################################################################
"""
    sif_simulation!(
                proj::ClumpingFactor2020{FT},
                year::FT
    ) where {FT<:AbstractFloat}

Simulate the SIF for Niwot Ridge, given
- `proj` [`ClumpingFactor2020`](@ref) type project control
- `year` Which year of data to simulate
- `ci` Clumping index
"""
function sif_simulation!(
            proj::ClumpingFactor2020{FT},
            year::Int,
            ci::FT
) where {FT<:AbstractFloat}
    @info tinfo("Simulating SIF for year $(year)...");

    # read the data from nc file to a CSV file, UTC correction made
    # please contact the leading author for the data
    _fold = "/home/wyujie/RAID/Data/ClumpingIndex";
    _file = "$(_fold)/sif_mip_climate_$(year)_ver.2018.11.28.nc";
    _df = DataFrame();
    _df[!,"Month" ]  = read_nc(_file, "month");
    _df[!,"Day"   ]  = read_nc(_file, "day");
    _df[!,"Hour"  ]  = read_nc(_file, "hour") .- 6;
    _df[!,"SW_in" ]  = read_nc(_file, "sw");
    _df[!,"SIF740"] .= NaN;
    _df[!,"SIF757"] .= NaN;
    _df[!,"SIF771"] .= NaN;

    # 0.1 unpack data
    _node = create_spac(proj, ci);
    @unpack angles, can_opt, can_rad, canopy_rt, envirs,f_SL,  ga, in_rad,
            latitude, leaves_rt, n_canopy, photo_set, plant_hs, plant_ps,
            rt_con, rt_dim, soil_opt, stomata_model, wl_set = _node;
    @unpack lidf, nAzi, nIncl = canopy_rt;
    @unpack dWL, iPAR = wl_set;
    _in_rad_bak = deepcopy(in_rad);
    _nSL = nAzi * nIncl;
    in_rad_sw = numerical∫(_in_rad_bak.E_direct , dWL) / 1000 +
                numerical∫(_in_rad_bak.E_diffuse, dWL) / 1000;
    angles.vza = 0;
    angles.raa = 0;

    # iterate through the weather data
    for i in eachindex(_df.Day)
        # update PAR related information
        zenith = zenith_angle(latitude, FT(_df.Day[i]), FT(_df.Hour[i]), FT(0));
        zenith = min(88, zenith);
        angles.sza = zenith;
        in_rad.E_direct  .= _in_rad_bak.E_direct  .* _df.SW_in[i] ./ in_rad_sw;
        in_rad.E_diffuse .= _in_rad_bak.E_diffuse .* _df.SW_in[i] ./ in_rad_sw;
        canopy_geometry!(canopy_rt, angles, can_opt, rt_con);
        canopy_matrices!(leaves_rt, can_opt);
        short_wave!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, rt_con);
        canopy_fluxes!(canopy_rt, can_opt, can_rad, in_rad, soil_opt,
                       leaves_rt, wl_set, rt_con);
        SIF_fluxes!(leaves_rt, can_opt, can_rad, canopy_rt, soil_opt, wl_set,
                    rt_con, rt_dim);

        # update SIF from the table
        _df.SIF740[i] = SIF_740(can_rad, wl_set);
        _df.SIF757[i] = SIF_757(can_rad, wl_set);
        _df.SIF771[i] = SIF_771(can_rad, wl_set);
    end;

    return _df
end



function sif_simulation!(
            proj::ClumpingFactor2020{FT}
) where {FT<:AbstractFloat}
    _fold = "/home/wyujie/RAID/Data/ClumpingIndex";
    _dfcs = sif_simulation!.([proj], 1999:2018, FT(0.48));
    save_csv!(vcat(_dfcs...), "$(_fold)/simulations/SIF_CI.csv");
    _dfcs = sif_simulation!.([proj], 1999:2018, FT(1));
    save_csv!(vcat(_dfcs...), "$(_fold)/simulations/SIF_nCI.csv");

    return nothing
end
