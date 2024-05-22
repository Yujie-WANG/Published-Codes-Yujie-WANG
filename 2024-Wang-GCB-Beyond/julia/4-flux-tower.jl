#
# this script is meant to run CliMA Land at the site level to examine the improvements in surface albedo
#
using DataFrames: DataFrame
using LazyArtifacts

using Emerald.EmeraldIO.Text: read_csv, save_csv!
using Emerald.EmeraldPhysics.EarthGeometry: solar_zenith_angle
using Emerald.EmeraldLand.Namespace: BulkSPAC, SPACConfiguration
using Emerald.EmeraldLand.SPAC: initialize_spac!, prescribe_soil!, prescribe_traits!, soil_plant_air_continuum!
using Emerald.EmeraldUtility.Time: parse_timestamp

FT = Float64;


# 1. read the data
flux_tower_df = read_csv("/home/wyujie/DATASERVER/field/AmeriFlux/original/AMF_US-MOz_BASE_HH_9-5.csv"; skiprows = 2);

# 2. define SPAC to work with
config = SPACConfiguration(FT);
spac = BulkSPAC(config);
prescribe_traits!(config, spac; cab = 57.23, ci = 0.69, lai = 4, sai = 1);
initialize_spac!(config, spac);
mask_uv = config.SPECTRA.Λ .< 400;
mask_par_nir = config.SPECTRA.Λ .>= 400;
rad_all = (config.SPECTRA.SOLAR_RAD[:,1] .+ config.SPECTRA.SOLAR_RAD[:,2])' * config.SPECTRA.ΔΛ * 1e-3;
rad_par_nir = (config.SPECTRA.SOLAR_RAD[mask_par_nir,1] .+ config.SPECTRA.SOLAR_RAD[mask_par_nir,2])' * config.SPECTRA.ΔΛ[mask_par_nir] * 1e-3;
spac_uv = deepcopy(spac);
spac_nouv = deepcopy(spac);

# 3. run the simulation for the entire July
mask_july = 201907000000 .< flux_tower_df.TIMESTAMP_START .< 201908000000;
sub_df = DataFrame();
sub_df[!,"TIME"] = flux_tower_df.TIMESTAMP_START[mask_july] .+ 15;
sub_df[!,"SW_IN"] = flux_tower_df.SW_IN_1_1_1[mask_july];
sub_df[!,"SWC"] = flux_tower_df.SWC_1_1_1[mask_july];
sub_df[!,"SW_OUT"] = flux_tower_df.SW_OUT_1_1_1[mask_july];
sub_df[!,"SW_OUT_UV"] .= 0.0;
sub_df[!,"SW_OUT_NOUV"] .= 0.0;
for dfr in eachrow(sub_df)
    spac_uv.meteo.rad_sw.e_dir .= dfr.SW_IN / rad_all .* config.SPECTRA.SOLAR_RAD[:,1];
    spac_uv.meteo.rad_sw.e_dif .= dfr.SW_IN / rad_all .* config.SPECTRA.SOLAR_RAD[:,2];
    spac_nouv.meteo.rad_sw.e_dir[mask_uv] .= 0;
    spac_nouv.meteo.rad_sw.e_dif[mask_uv] .= 0;
    spac_nouv.meteo.rad_sw.e_dir[mask_par_nir] .= dfr.SW_IN / rad_par_nir .* config.SPECTRA.SOLAR_RAD[mask_par_nir,1];
    spac_nouv.meteo.rad_sw.e_dif[mask_par_nir] .= dfr.SW_IN / rad_par_nir .* config.SPECTRA.SOLAR_RAD[mask_par_nir,2];
    # run the model for 0 seconds so as to compute the shortwave radiation only
    fdoy = parse_timestamp(dfr.TIME; in_format = "YYYYMMDDhhmm", out_format = "FDOY");
    sza = solar_zenith_angle(38.74, fdoy);
    spac_uv.canopy.sun_geometry.state.sza = (dfr.SW_IN > 1) ? min(sza, 88.999) : sza;
    spac_nouv.canopy.sun_geometry.state.sza = (dfr.SW_IN > 1) ? min(sza, 88.999) : sza;
    swc = dfr.SWC / 100;
    swcs = Tuple([swc for i in 1:length(spac.soils)]);
    prescribe_soil!(spac_uv; swcs = swcs);
    soil_plant_air_continuum!(config, spac_uv, 0);
    soil_plant_air_continuum!(config, spac_nouv, 0);
    dfr.SW_OUT_UV = spac_uv.canopy.sun_geometry.auxil.e_difꜛ[:,1]' * config.SPECTRA.ΔΛ * 1e-3;
    dfr.SW_OUT_NOUV = spac_nouv.canopy.sun_geometry.auxil.e_difꜛ[:,1]' * config.SPECTRA.ΔΛ * 1e-3;
end;

# 4. save the simulation results
save_csv!(sub_df, "../output/4_flux_tower.csv");
