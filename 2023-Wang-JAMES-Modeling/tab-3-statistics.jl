using DataFrames: DataFrame

using GriddingMachine.Indexer: lat_ind, lon_ind

using JuliaUtilities.EmeraldRegression: linear_regress
using JuliaUtilities.NetcdfIO: read_nc
using JuliaUtilities.TextIO: save_csv!


CLIMA_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
CLIMA_TAG = "a6_gm1_wd1";
FIG_FOLDER = "";


function generate_df!()
    _nc_file = "$(CLIMA_FOLDER)/$(CLIMA_TAG)/$(CLIMA_TAG)_2019_1X_1D.c3.rpar.nc";

    # create data frame
    _df = DataFrame();
    _df[!,"ID"] = ["a", "b", "c", "d", "e", "f", "g", "h"];
    _df[!,"Site"] = ["US-NR1", "BR-Npw", "GH-Ank", "DE-Hai", "ZA-Kru", "RU-Skp", "CN-Qia", "AU-Tum"];
    _df[!,"GPP-SIF683" ] .= NaN;
    _df[!,"GPP-SIF740" ] .= NaN;
    _df[!,"GPP-SIF757" ] .= NaN;
    _df[!,"GPP-SIF771" ] .= NaN;
    _df[!,"GPP-NDVI"   ] .= NaN;
    _df[!,"GPP-EVI"    ] .= NaN;
    _df[!,"GPP-NIRv"   ] .= NaN;
    _df[!,"SIF683-NIRv"] .= NaN;
    _df[!,"SIF740-NIRv"] .= NaN;
    _df[!,"SIF757-NIRv"] .= NaN;
    _df[!,"SIF771-NIRv"] .= NaN;

    # iterate through the
    _lats = [  40.0329, -16.4980,  5.2685, 51.0792, -25.0197,  62.2550,  26.7414, -35.6566];
    _lons = [-105.5464, -56.4120, -2.6942, 10.4522,  31.4969, 129.1680, 115.0581, 148.1517];
    _lati = lat_ind.(_lats; res=1.0);
    _loni = lon_ind.(_lons; res=1.0);
    for _i in eachindex(_lats)
        _gpp    = read_nc(_nc_file, "mGPP"   , _loni[_i], _lati[_i]);
        _sif683 = read_nc(_nc_file, "mSIF683", _loni[_i], _lati[_i]);
        _sif740 = read_nc(_nc_file, "mSIF740", _loni[_i], _lati[_i]);
        _sif757 = read_nc(_nc_file, "mSIF757", _loni[_i], _lati[_i]);
        _sif771 = read_nc(_nc_file, "mSIF771", _loni[_i], _lati[_i]);
        _ndvi   = read_nc(_nc_file, "mNDVI"  , _loni[_i], _lati[_i]);
        _evi    = read_nc(_nc_file, "mEVI"   , _loni[_i], _lati[_i]);
        _nirv   = read_nc(_nc_file, "mNIRvI" , _loni[_i], _lati[_i]);

        # GPP ~ SIF
        _lr1  = linear_regress((_sif683[:],1), _gpp[:]   );
        _lr2  = linear_regress((_sif740[:],1), _gpp[:]   );
        _lr3  = linear_regress((_sif757[:],1), _gpp[:]   );
        _lr4  = linear_regress((_sif771[:],1), _gpp[:]   );
        _lr5  = linear_regress((_ndvi[:]  ,1), _gpp[:]   );
        _lr6  = linear_regress((_evi[:]   ,1), _gpp[:]   );
        _lr7  = linear_regress((_nirv[:]  ,1), _gpp[:]   );
        _lr8  = linear_regress((_nirv[:]  ,1), _sif683[:]);
        _lr9  = linear_regress((_nirv[:]  ,1), _sif740[:]);
        _lr10 = linear_regress((_nirv[:]  ,1), _sif757[:]);
        _lr11 = linear_regress((_nirv[:]  ,1), _sif771[:]);
        _df[_i,3 ] = round(_lr1.r2 ; digits=3);
        _df[_i,4 ] = round(_lr2.r2 ; digits=3);
        _df[_i,5 ] = round(_lr3.r2 ; digits=3);
        _df[_i,6 ] = round(_lr4.r2 ; digits=3);
        _df[_i,7 ] = round(_lr5.r2 ; digits=3);
        _df[_i,8 ] = round(_lr6.r2 ; digits=3);
        _df[_i,9 ] = round(_lr7.r2 ; digits=3);
        _df[_i,10] = round(_lr8.r2 ; digits=3);
        _df[_i,11] = round(_lr9.r2 ; digits=3);
        _df[_i,12] = round(_lr10.r2; digits=3);
        _df[_i,13] = round(_lr11.r2; digits=3);
    end;

    # save the data frame
    save_csv!(_df, "$(FIG_FOLDER)/2021_clima_land_tab_3_statistics.csv");

    return nothing
end;
