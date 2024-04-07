###############################################################################
#
# Tower sites
#
###############################################################################
""" abstract type AbstractFluxTower """
abstract type AbstractFluxTower end




"""
    struct Tower_US_MOz end

Flux tower at Ozark, Missouri, USA. This is a deciduous forest (C3 angiosperm).
    Website: https://ameriflux.lbl.gov/sites/siteinfo/US-MOz
"""
struct Tower_US_MOz end




"""
    struct Tower_US_NR1

Flux tower at Niwot Ridge, Colorado, USA, This is a subalphine evergreen forest
    (C3 gymnosperm). Website: https://ameriflux.lbl.gov/sites/siteinfo/US-NR1
"""
struct Tower_US_NR1 end








###############################################################################
#
# Function to reprocess CSV files to meet our project requirements
#
###############################################################################
"""
Reprocess FLUXNET2015 CSV files to meet project requirements.
"""
function reprocess_FLUXNET! end




"""
This method reprocesses the data for `US_MOz` site:

    reprocess_FLUXNET!(site::Tower_US_MOz)

Reprocess FLUXNET2015 data to per year CSV files, given
- `site` [`Tower_US_MOz`](@ref) type site indentifier
"""
reprocess_FLUXNET!(site::Tower_US_MOz) =
(
    # location of the files
    _file_loc = "/home/wyujie/RAID/Data/FLUXNET2015/";

    # columns needed for Ozark:
    #     TIMESTAMP_START   => TIME
    #     PPFD_IN_1_1_1     => PPFD
    #     LW_OUT_1_1_1      => LW_OUT
    #     TA_1_1_1          => T_AIR
    #     WS_1_1_1          => WIND
    #     TS_1_1_1          => T_SOIL
    #     SWC_1_1_1         => SWC
    #     PA_1_1_1          => P_ATM
    #     FC_1_1_1          => F_CO2
    #     LE_1_1_1          => LE_H2O
    #     RH_1_1_1          -> RH,VPD
    # need to filter and rename the columns
    _data = read_csv("$(_file_loc)US_MOz.csv"; skiprows=2);

    # mask the unrealistic numbers to NaN
    for col in 3:length(names(_data))
        _mask = _data[:,col] .< -9000;
        _data[_mask,col] .= NaN;
    end;

    # copy data to new DataFrame
    _eata = DataFrame();
    _eata.TIME   = _data.TIMESTAMP_START;
    _eata.F_CO2  = _data.FC_1_1_1;
    _eata.LE_H2O = _data.LE_1_1_1;
    _eata.LW_OUT = _data.LW_OUT_1_1_1;
    _eata.P_ATM  = _data.PA_1_1_1;
    _eata.PPFD   = _data.PPFD_IN_1_1_1;
    _eata.SWC    = _data.SWC_1_1_1;
    _eata.T_AIR  = _data.TA_1_1_1;
    _eata.T_SOIL = _data.TS_1_1_1;
    _eata.RH     = _data.RH_1_1_1;
    _eata.WIND   = _data.WS_1_1_1;

    # mask the unrealistic numbers to NaN
    for _col in 2:length(names(_eata))
        _mask = _eata[:,_col] .< -9000;
        _eata[_mask,_col] .= NaN;
    end;
    _eata.VPD = saturation_vapor_pressure.(_eata.T_AIR .+ 273.15) ./100 .*
                (100 .- _data.RH_1_1_1) ./ 100;

    # unselect any row with a NaN
    _mask = isnan.(_eata.F_CO2 ) .| isnan.(_eata.LE_H2O) .|
            isnan.(_eata.LW_OUT) .| isnan.(_eata.P_ATM ) .|
            isnan.(_eata.PPFD  ) .| isnan.(_eata.RH    ) .|
            isnan.(_eata.SWC   ) .| isnan.(_eata.T_AIR ) .|
            isnan.(_eata.T_SOIL) .| isnan.(_eata.VPD   ) .|
            isnan.(_eata.WIND );
    _fata = _eata[.!_mask,:];

    # save the data per year
    for _year in 2006:2020
        _mask = (_fata.TIME .> _year*1e8) .* (_fata.TIME .< (_year+1)*1e8);
        if sum(_mask) > 10000
            save_csv!("$(_file_loc)Ozark_$(_year).csv", _fata[_mask, :]);
        end;
    end;

    return nothing
)




"""
This method reprocesses the data for `US_NR1` site:

    reprocess_FLUXNET!(site::Tower_US_MOz)

Reprocess FLUXNET2015 data to per year CSV files, given
- `site` [`Tower_US_NR1`](@ref) type site indentifier
"""
reprocess_FLUXNET!(site::Tower_US_NR1) =
(
    # location of the files
    _file_loc = "/home/wyujie/RAID/Data/FLUXNET2015/";

    # columns needed for Niwot Ridge:
    #     FC_PI_F_1_1_1         | carbon flux               => F_CO2
    #     LE_PI_F_1_1_1         | latent heat flux          => F_H2O
    #     LW_OUT_PI_F_1_1_1     | leaf temperature          => LW_OUT
    #     PA_PI_F_1_1_1         | atmospheric pressure      => P_ATM
    #     PPFD_IN_PI_F_1_1_1    | correct in_rad            => PPFD
    #     SWC_PI_F_1_1_1        | for p_soil                => SWC
    #     TA_PI_F_1_1_1         |--------------------------|
    #     TA_PI_F_1_2_1         |--------------------------|
    #     TA_PI_F_1_3_1         | mean air temperature-----|=> T_AIR
    #     TIMESTAMP_START       | time                      => TIME_START
    #     TS_PI_F_1_1_1         | soil temperature          => T_SOIL
    #     VPD_PI_F_1_1_1        |--------------------------|
    #     VPD_PI_F_1_2_1        |--------------------------|
    #     VPD_PI_F_1_3_1        | vapor pressure deficit---|=> VPD
    #     WS_PI_F_1_1_1         | wind speed                => WIND
    # need to filter and rename the columns
    _data = read_csv("$(_file_loc)US_NR1.csv"; skiprows=2);

    # mask the unrealistic numbers to NaN
    for _col in 3:length(names(_data))
        _mask = _data[:,_col] .< -9000;
        _data[_mask,_col] .= NaN;
    end;

    # copy data to new DataFrame
    _eata = DataFrame();
    _eata.TIME   = _data.TIMESTAMP_START;
    _eata.F_CO2  = _data.FC_PI_F_1_1_1;
    _eata.LE_H2O = _data.LE_PI_F_1_1_1;
    _eata.LW_OUT = _data.LW_OUT_PI_F_1_1_1;
    _eata.P_ATM  = _data.PA_PI_F_1_1_1;
    _eata.PPFD   = _data.PPFD_IN_PI_F_1_1_1;
    _eata.SWC    = _data.SWC_PI_F_1_1_1;
    _eata.T_AIR  = [nanmean([_data.TA_PI_F_1_1_1[_i],
                             _data.TA_PI_F_1_2_1[_i],
                             _data.TA_PI_F_1_3_1[_i]])
                    for _i in eachindex(_data.TIMESTAMP_START)];
    _eata.T_SOIL = _data.TS_PI_F_1_1_1;
    _eata.VPD    = [nanmean([_data.VPD_PI_F_1_1_1[_i],
                             _data.VPD_PI_F_1_2_1[_i],
                             _data.VPD_PI_F_1_3_1[_i]])
                    for _i in eachindex(_data.TIMESTAMP_START)];
    _eata.WIND   = _data.WS_PI_F_1_1_1;
    _eata.RH     = 100 .- 100 .* _eata.VPD .* 100 ./
                          saturation_vapor_pressure.(_eata.T_AIR .+ 273.15);

    # mask the unrealistic numbers to NaN
    for _col in 2:length(names(_eata))
        _mask = _eata[:,_col] .< -9000;
        _eata[_mask,_col] .= NaN;
    end;

    # unselect any row with a NaN
    _mask = isnan.(_eata.F_CO2 ) .| isnan.(_eata.LE_H2O) .|
            isnan.(_eata.LW_OUT) .| isnan.(_eata.P_ATM ) .|
            isnan.(_eata.PPFD  ) .| isnan.(_eata.RH    ) .|
            isnan.(_eata.SWC   ) .| isnan.(_eata.T_AIR ) .|
            isnan.(_eata.T_SOIL) .| isnan.(_eata.VPD   ) .|
            isnan.(_eata.WIND  );
    _fata = _eata[.!_mask,:];

    # save the data per year
    for _year in 2006:2020
        _mask = (_fata.TIME .> _year*1e8) .* (_fata.TIME .< (_year+1)*1e8);
        if sum(_mask) > 10000
            save_csv!("$(_file_loc)NiwotRidge_$(_year).csv", _fata[_mask, :]);
        end;
    end;

    return nothing
)
