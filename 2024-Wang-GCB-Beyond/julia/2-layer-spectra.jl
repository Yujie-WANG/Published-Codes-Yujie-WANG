#
# this script is used to output the radiation spectra at different layers of the canopy
#
using DataFrames: DataFrame

using Emerald.EmeraldIO.Text: save_csv!
using Emerald.EmeraldLand.Namespace: BulkSPAC, LAND_2021_1NM, SPACConfiguration
using Emerald.EmeraldLand.SPAC: initialize_spac!, soil_plant_air_continuum!

FT = Float64;


# 1 run the model at two contrasting configurations
config = SPACConfiguration(FT, dataset = LAND_2021_1NM);
spac = BulkSPAC(config);
initialize_spac!(config, spac);
soil_plant_air_continuum!(config, spac, FT(1));

# 2 get the radiation spectra at different layers of the canopy
df = DataFrame(WL = config.SPECTRA.Λ);
for irt in eachindex(spac.plant.leaves)
    rad = spac.canopy.sun_geometry.auxil.e_dirꜜ[:,irt] .+ spac.canopy.sun_geometry.auxil.e_difꜜ[:,irt] .+ spac.canopy.sun_geometry.auxil.e_difꜛ[:,irt];
    df[!,"SPECTRUM_$(irt)"] = rad;
end;

# 3 save the data
save_csv!("../output/2_radiation_spectra.csv", df);
