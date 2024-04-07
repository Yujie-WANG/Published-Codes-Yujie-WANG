#
# this script is meant for generating data to train empirical models for the f_ppar
#
using DataFrames: DataFrame
using ProgressMeter: @showprogress

using Emerald.EmeraldIO.Text: save_csv!
using Emerald.EmeraldLand.Namespace: BulkSPAC, LAND_2021_1NM, SPACConfiguration
using Emerald.EmeraldLand.SPAC: initialize!, prescribe_traits!, soil_plant_air_continuum!
using Emerald.EmeraldMath.Stats: nanmean
using Emerald.EmeraldMath.Regression: linear_regress

FT = Float64;


# 1. generate random canopy f_ppar
config = SPACConfiguration{FT}(DATASET = LAND_2021_1NM, ENABLE_REF = false, ENABLE_SIF = false);
spac = BulkSPAC(config);
initialize!(config, spac);
soil_plant_air_continuum!(config, spac, FT(0));
df = DataFrame(LAI = Float64[], CHL = Float64[], Z = Float64[], F_PPAR_SL = Float64[], F_PPAR_SH = Float64[]);
@showprogress for i in 1:1000
    lai = rand(0.5:0.1:6.0);
    chl = rand(5:1:80);
    prescribe_traits!(config, spac, lai = lai, cab = chl, car = chl / 7);
    soil_plant_air_continuum!(config, spac, FT(0));
    for ilf in eachindex(spac.plant.leaves)
        apar_sl = spac.plant.leaves[ilf].flux.auxil.apar_sunlit;
        apar_sh = spac.plant.leaves[ilf].flux.auxil.apar_shaded;
        ppar_sl = spac.plant.leaves[ilf].flux.auxil.ppar_sunlit;
        ppar_sh = spac.plant.leaves[ilf].flux.auxil.ppar_shaded;
        push!(df, (lai, chl, ilf / length(spac.plant.leaves), nanmean(ppar_sl ./ apar_sl), ppar_sh / apar_sh));
    end;
end;

# 2. run the regression
lr_sl = linear_regress((df.Z, df.LAI, log.(log.(df.CHL)), 1), df.F_PPAR_SL);
lr_sh = linear_regress((df.Z, df.LAI, log.(log.(df.CHL)), 1), df.F_PPAR_SH);
df[!,"F_PPAR_SL_PRED"] .= lr_sl.XY.predY;
df[!,"F_PPAR_SH_PRED"] .= lr_sh.XY.predY;
#=
lr_sl.LM:
─────────────────────────────────────────────────────────────────────────
          Coef.   Std. Error        t  Pr(>|t|)    Lower 95%    Upper 95%
─────────────────────────────────────────────────────────────────────────
x1  -0.0199109   0.00015179   -131.17    <1e-99  -0.0202085   -0.0196134
x2   0.00123618  2.73773e-5     45.15    <1e-99   0.00118252   0.00128985
x3   0.210514    0.000187843  1120.69    <1e-99   0.210146     0.210882
x4   0.510661    0.000266471  1916.38    <1e-99   0.510138     0.511183
─────────────────────────────────────────────────────────────────────────
lr_sh.LM:
─────────────────────────────────────────────────────────────────────────
          Coef.   Std. Error        t  Pr(>|t|)    Lower 95%    Upper 95%
─────────────────────────────────────────────────────────────────────────
x1  -0.0550177   0.000257525  -213.64    <1e-99  -0.0555225   -0.0545129
x2  -0.00309678  4.64481e-5    -66.67    <1e-99  -0.00318783  -0.00300574
x3   0.175522    0.000318692   550.76    <1e-99   0.174898     0.176147
x4   0.549334    0.000452092  1215.09    <1e-99   0.548448     0.55022
─────────────────────────────────────────────────────────────────────────
=#

# 3. save the results
save_csv!("../output/6_empirical.csv", df);

# 4. the case with combined sunlit and shaded leaves
new_df = DataFrame(LAI = [df.LAI; df.LAI], CHL = [df.CHL; df.CHL], Z = [df.Z; df.Z], F_PPAR = [df.F_PPAR_SL; df.F_PPAR_SH]);
new_lr = linear_regress((new_df.Z, new_df.LAI, log.(log.(new_df.CHL)), 1), new_df.F_PPAR);
#=
new_lr.LM:
    f_ppar = -0.0374643 * z + -0.0009303 * lai + 0.193018 * log(log(chl)) + 0.529997
────────────────────────────────────────────────────────────────────────
         Coef.   Std. Error       t  Pr(>|t|)    Lower 95%     Upper 95%
────────────────────────────────────────────────────────────────────────
x1  -0.0374643  0.000481077  -77.88    <1e-99  -0.0384073   -0.0365214
x2  -0.0009303  8.67685e-5   -10.72    <1e-26  -0.00110037  -0.000760228
x3   0.193018   0.000595341  324.21    <1e-99   0.191851     0.194185
x4   0.529997   0.000844542  627.56    <1e-99   0.528342     0.531653
────────────────────────────────────────────────────────────────────────
=#
