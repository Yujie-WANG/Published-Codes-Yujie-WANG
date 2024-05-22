#
# this script is meant to run CliMA Land at the site level to examine the improvements in surface albedo
#
using Emerald.EmeraldIO.Text: read_csv
using Emerald.EmeraldMath.Regression: linear_regress

FT = Float64;


# 1. read the simulation results
df = read_csv("../output/4_flux_tower.csv");
mask = (df.SW_OUT_UV .> 0) .&& (df.SW_OUT_NOUV .> 0);

lr3 = linear_regress((df.SW_OUT[mask],1), df.SW_OUT_UV[mask]);
lr4 = linear_regress((df.SW_OUT[mask],1), df.SW_OUT_NOUV[mask]);
@show lr3.LM;
@show lr4.LM;

lr3_0 = linear_regress((df.SW_OUT[mask],), df.SW_OUT_UV[mask]);
lr4_0 = linear_regress((df.SW_OUT[mask],), df.SW_OUT_NOUV[mask]);
@show lr3_0.LM;
@show lr4_0.LM;
