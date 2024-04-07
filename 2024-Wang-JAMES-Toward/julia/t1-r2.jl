#
# this script is meant for computing how much the bias in GPP can be explained by the tested parameters
#
using DataFrames: DataFrame

using Emerald.EmeraldIO.Text: save_csv!
using Emerald.EmeraldMath.Regression: linear_regress
using Emerald.EmeraldMath.Stats: nanmean, nanstd
using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: read_LUT
using NetcdfIO: read_nc


# 1. read the data from the netcdf file
ncfile = "../output/4_gpp_sif.nc";
delta_gpp = read_nc(ncfile, "GPP_BB_KDST") .- read_nc(ncfile, "GPP_HS_NVER");
delta_sif = read_nc(ncfile, "SIF_BB_KDST") .- read_nc(ncfile, "SIF_HS_KDTD");

# 2. read the global maps
chl_7d = regrid(read_LUT(query_collection("CHL_2X_7D_V1"))[1], 1);
cli_1m = regrid(read_LUT(query_collection("CI_2X_1M_V3"))[1], 1);
lai_8d = regrid(read_LUT(query_collection("LAI_MODIS_2X_8D_2019_V1"))[1], 1);
sla_1y = regrid(read_LUT(query_collection("SLA_2X_1Y_V1"))[1], 1);
chl_1y = ones(360,180) .* NaN;
cli_1y = ones(360,180) .* NaN;
lai_1y = ones(360,180) .* NaN;
lma_1y = 1 ./ sla_1y;
for ilon in 1:360, ilat in 1:180
    chl_1y[ilon,ilat] = nanmean(chl_7d[ilon,ilat,:]);
    cli_1y[ilon,ilat] = nanmean(cli_1m[ilon,ilat,:]);
    lai_1y[ilon,ilat] = nanmean(lai_8d[ilon,ilat,:]);
end;

# 3. create a dataframe to store the results
df = DataFrame();
lbs = ["Y ~ CHL",
       "Y ~ CI",
       "Y ~ LAI",
       "Y ~ LMA",
       "Y ~ LAI + CHL",
       "Y ~ LAI + CI",
       "Y ~ LAI + LMA",
       "Y ~ CHL + CI",
       "Y ~ CHL + LMA",
       "Y ~ CI + LMA",
       "Y ~ LAI + CHL + CI",
       "Y ~ LAI + CHL + LMA",
       "Y ~ LAI + CI + LMA",
       "Y ~ CHL + CI + LMA",
       "Y ~ LAI + CHL + CI + LMA"];
df[!, "XY"   ]  = lbs;
df[!, "R²GPP"] .= NaN;
df[!, "R²SIF"] .= NaN;

# 4. run the regression for multiple variables
ys = [delta_gpp, delta_sif];
for i in 1:2
    lr = linear_regress((chl_1y[:],1), ys[i][:]);
    label = i==1 ? "R²GPP" : "R²SIF";
    df[1, label] = lr.R²;

    lr = linear_regress((cli_1y[:],1), ys[i][:]);
    df[2, label] = lr.R²;

    lr = linear_regress((lai_1y[:],1), ys[i][:]);
    df[3, label] = lr.R²;

    lr = linear_regress((lma_1y[:],1), ys[i][:]);
    df[4, label] = lr.R²;

    lr = linear_regress((lai_1y[:],chl_1y[:],1), ys[i][:]);
    df[5, label] = lr.R²;

    lr = linear_regress((lai_1y[:],cli_1y[:],1), ys[i][:]);
    df[6, label] = lr.R²;

    lr = linear_regress((lai_1y[:],lma_1y[:],1), ys[i][:]);
    df[7, label] = lr.R²;

    lr = linear_regress((chl_1y[:],cli_1y[:],1), ys[i][:]);
    df[8, label] = lr.R²;

    lr = linear_regress((chl_1y[:],lma_1y[:],1), ys[i][:]);
    df[9, label] = lr.R²;

    lr = linear_regress((cli_1y[:],lma_1y[:],1), ys[i][:]);
    df[10, label] = lr.R²;

    lr = linear_regress((lai_1y[:],chl_1y[:],cli_1y[:],1), ys[i][:]);
    df[11, label] = lr.R²;

    lr = linear_regress((lai_1y[:],chl_1y[:],lma_1y[:],1), ys[i][:]);
    df[12, label] = lr.R²;

    lr = linear_regress((lai_1y[:],cli_1y[:],lma_1y[:],1), ys[i][:]);
    df[13, label] = lr.R²;

    lr = linear_regress((chl_1y[:],cli_1y[:],lma_1y[:],1), ys[i][:]);
    df[14, label] = lr.R²;

    lr = linear_regress((lai_1y[:],chl_1y[:],cli_1y[:],lma_1y[:],1), ys[i][:]);
    df[15, label] = lr.R²;
end;

# 5. save the results
save_csv!(df, "../output/t1_r2.csv");
