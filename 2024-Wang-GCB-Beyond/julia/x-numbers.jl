#
# This script is meant to calculate the numbers to use in the main text
#
using Emerald.EmeraldLand.Namespace: BulkSPAC, LAND_2021_1NM, SPACConfiguration
using Emerald.EmeraldLand.SPAC: initialize_spac!
using Emerald.EmeraldPhysics.Optics: photon

FT = Float64;


# 1. create the configuration
config = SPACConfiguration(FT; dataset = LAND_2021_1NM);
srad_energy = config.SPECTRA.SOLAR_RAD[:,1] + config.SPECTRA.SOLAR_RAD[:,2];
srad_photon = photon.(config.SPECTRA.Λ, srad_energy);

# 2. calculate the number of photons in UV, PAR, and FR
mask_uv = (config.SPECTRA.Λ .< 400);
mask_par = (400 .<= config.SPECTRA.Λ .<= 700);
mask_fr = (700 .< config.SPECTRA.Λ .<=750);

@show srad_energy[mask_uv]' * config.SPECTRA.ΔΛ[mask_uv] * 1e-3;
@show srad_energy[mask_par]' * config.SPECTRA.ΔΛ[mask_par] * 1e-3;
@show srad_energy[mask_fr]' * config.SPECTRA.ΔΛ[mask_fr] * 1e-3;

@show srad_photon[mask_uv]' * config.SPECTRA.ΔΛ[mask_uv] * 1e-3 * 1e6;
@show srad_photon[mask_par]' * config.SPECTRA.ΔΛ[mask_par] * 1e-3 * 1e6;
@show srad_photon[mask_fr]' * config.SPECTRA.ΔΛ[mask_fr] * 1e-3 * 1e6;

# 3. some information about the SPAC
spac = BulkSPAC(config);
initialize_spac!(config, spac);
