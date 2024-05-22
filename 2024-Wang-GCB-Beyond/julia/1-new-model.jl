#
# This script is meant to generate the model output to illustrate that the new photosynthesis model is working as expected.
#
using DataFrames: DataFrame

using Emerald.EmeraldIO.Text: save_csv!
using Emerald.EmeraldLand.Namespace: LAND_2021_1NM, SPACConfiguration

FT = Float64;


# 1. create the configuration and spac structs for the simulation
config = SPACConfiguration(FT; dataset = LAND_2021_1NM);

# 2. save the natural radiation and absorption spectra of the key contents
save_csv!("../output/1_spectra.csv",
          DataFrame(WL        = config.SPECTRA.Λ,
                    E_DIR     = config.SPECTRA.SOLAR_RAD[:,1],
                    E_DIF     = config.SPECTRA.SOLAR_RAD[:,2],
                    REF_WATER = config.SPECTRA.K_H₂O,
                    REF_CHL   = config.SPECTRA.K_CAB,
                    REF_CARV  = config.SPECTRA.K_CAR_V,
                    REF_CARZ  = config.SPECTRA.K_CAR_Z,
                    REF_LMA   = config.SPECTRA.K_LMA));
