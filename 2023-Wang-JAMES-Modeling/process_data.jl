using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: read_LUT

using JuliaUtilities.EmeraldMath: nanmean


# load land mask
begin
    LAI_2X,_        = read_LUT(query_collection("LAI_MODIS_2X_1M_2019_V1"));
    LAI_1X          = regrid(LAI_2X, 1);
    LAIMASK_1X      = ones(Bool,360,180);
    LM_4X,_         = read_LUT(query_collection("LM_4X_1Y_V1"));
    LM_1X           = regrid(LM_4X, 1);
    MASK_1L         = (LM_1X .<= 0.1);
    LANDMASK_1M     = ones(Bool,360,180,12);
    LANDMASK_8D     = ones(Bool,360,180,46);
    LANDMASK_1M_LAI = ones(Bool,360,180,12);
    LANDMASK_8D_LAI = ones(Bool,360,180,46);
    for _lon in 1:360, _lat in 1:180
        LAIMASK_1X[_lon,_lat] = (maximum(LAI_1X[_lon,_lat,:]) < 2);
    end;
    for _l in 1:12
        LANDMASK_1M[:,:,_l] .= MASK_1L;
        LANDMASK_1M_LAI[:,:,_l] .= MASK_1L .|| LAIMASK_1X;
    end;
    for _l in 1:46
        LANDMASK_8D[:,:,_l] .= MASK_1L;
        LANDMASK_8D_LAI[:,:,_l] .= MASK_1L .|| LAIMASK_1X;
    end;
end;


# function to process the data to plot
function processed_data(ref_1x, clm_1x; laimask=false)
    # filter site without land to NaN
    if laimask
        if size(ref_1x,3) == 12
            ref_1x[LANDMASK_1M_LAI] .= NaN;
            clm_1x[LANDMASK_1M_LAI] .= NaN;
        elseif size(ref_1x,3) == 46
            ref_1x[LANDMASK_8D_LAI] .= NaN;
            clm_1x[LANDMASK_8D_LAI] .= NaN;
        end;
    else
        if size(ref_1x,3) == 12
            ref_1x[LANDMASK_1M] .= NaN;
            clm_1x[LANDMASK_1M] .= NaN;
        elseif size(ref_1x,3) == 46
            ref_1x[LANDMASK_8D] .= NaN;
            clm_1x[LANDMASK_8D] .= NaN;
        end;
    end;

    # calculate the annual mean
    _ref_mean = ones(360,180);
    _clm_mean = ones(360,180);
    for _ilon in 1:360, _ilat in 1:180
        _mask = .!(isnan.(ref_1x[_ilon,_ilat,:]) .|| isnan.(clm_1x[_ilon,_ilat,:]));
        _ref_mean[_ilon,_ilat] = nanmean(ref_1x[_ilon,_ilat,_mask]);
        _clm_mean[_ilon,_ilat] = nanmean(clm_1x[_ilon,_ilat,_mask]);
    end;

    # add extra filter to land coverage
    _mask = (LM_1X .< 0.8);
    _ref_mean[_mask] .= NaN;
    _clm_mean[_mask] .= NaN;

    # filter the data
    _mask = (_clm_mean .> 0);
    _data_1 = _ref_mean[_mask];
    _data_2 = _clm_mean[_mask];

    return _data_1, _data_2
end;


# function to process the data to plot (for global maps)
function processed_data_2d(ref_1x, clm_1x)
    # filter site without land to NaN
    if size(ref_1x,3) == 12
        ref_1x[LANDMASK_1M] .= NaN;
        clm_1x[LANDMASK_1M] .= NaN;
    elseif size(ref_1x,3) == 46
        ref_1x[LANDMASK_8D] .= NaN;
        clm_1x[LANDMASK_8D] .= NaN;
    end;

    # calculate the annual mean
    _ref_mean = ones(360,180);
    _clm_mean = ones(360,180);
    for _ilon in 1:360, _ilat in 1:180
        _mask = .!(isnan.(ref_1x[_ilon,_ilat,:]) .|| isnan.(clm_1x[_ilon,_ilat,:]));
        _ref_mean[_ilon,_ilat] = nanmean(ref_1x[_ilon,_ilat,_mask]);
        _clm_mean[_ilon,_ilat] = nanmean(clm_1x[_ilon,_ilat,_mask]);
    end;

    # add extra filter to land coverage
    _mask = (LM_1X .< 0.8);
    _ref_mean[_mask] .= NaN;
    _clm_mean[_mask] .= NaN;

    return _ref_mean, _clm_mean
end;


# function to process the data to plot
function processed_data_3d(ref_1x, clm_1x)
    # filter site without land to NaN
    if size(ref_1x,3) == 12
        ref_1x[LANDMASK_1M] .= NaN;
        clm_1x[LANDMASK_1M] .= NaN;
    elseif size(ref_1x,3) == 46
        ref_1x[LANDMASK_8D] .= NaN;
        clm_1x[LANDMASK_8D] .= NaN;
    end;

    return ref_1x, clm_1x
end;
