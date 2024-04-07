#
# this script is supposed to generate the spectra and absorption features of the leaf to show the problem with broadband models
#
using DataFrames: DataFrame

using Emerald.EmeraldIO.Text: save_csv!
using Emerald.EmeraldLand.LeafOptics: leaf_spectra!
using Emerald.EmeraldLand.Namespace: LAND_2021_1NM, Leaf, SPACConfiguration, ShortwaveRadiation
using Emerald.EmeraldPhysics.Optics: photon

FT = Float64;
config = SPACConfiguration{FT}(DATASET = LAND_2021_1NM, ENABLE_SIF = false);
leaf = Leaf(config);

# 1. output the reference spectra within 400-700 nm
save_csv!("../output/1_spectra.csv",
          DataFrame(WL = config.SPECTRA.Λ, KCHL = config.SPECTRA.K_CAB, KCARV = config.SPECTRA.K_CAR_V, KCARZ = config.SPECTRA.K_CAR_Z, KH2O = config.SPECTRA.K_H₂O, KLMA = config.SPECTRA.K_LMA));

# 2. create the radiation and leaf objects
rad_nature = ShortwaveRadiation(config);
rad_zero = deepcopy(rad_nature);
rad_zero.e_dir .= 0;
rad_zero.e_dif .= 0;
rad_blue = deepcopy(rad_zero);
rad_blue.e_dir[448 .<= config.SPECTRA.Λ .<= 458] .= 5000;
rad_green = deepcopy(rad_zero);
rad_green.e_dir[530 .<= config.SPECTRA.Λ .<= 540] .= 5000;
rad_red = deepcopy(rad_zero);
rad_red.e_dir[655 .<= config.SPECTRA.Λ .<= 665] .= 5000;

# 3. define function to compute PAR, APAR, and PPAR from Leaf and ShortwaveRadiation
function leaf_fapar_fppar(config::SPACConfiguration{FT}, leaf::Leaf{FT}, rad::ShortwaveRadiation{FT}) where {FT}
    (; IΛ_PAR, ΔΛ_PAR, Λ) = config.SPECTRA;

    # convert the radiation to photon flux
    photon_dir = photon.(Λ, rad.e_dir) .* 1000;
    photon_dif = photon.(Λ, rad.e_dif) .* 1000;
    photon_all = photon_dir .+ photon_dif;

    # compute the absorbed PAR by the leaf
    photon_leaf = leaf.bio.auxil.α_leaf .* photon_all;

    # compute the PAR that goes to photosystems
    photon_ppar = leaf.bio.auxil.f_ppar .* photon_leaf;

    # compute the sum of par, apar, and ppar
    par = photon_all[IΛ_PAR]' * ΔΛ_PAR;
    apar = photon_leaf[IΛ_PAR]' * ΔΛ_PAR;
    ppar = photon_ppar[IΛ_PAR]' * ΔΛ_PAR;

    return apar / par, ppar / apar
end;

# 4. compute the PAR, APAR, and PPAR for the leaf with different chorophyll contents at different radiation conditions
df = DataFrame(CHL = collect(FT,1:80),
               F_APAR_NATURE = zeros(FT,80),
               F_APAR_BLUE = zeros(FT,80),
               F_APAR_GREEN = zeros(FT,80),
               F_APAR_RED = zeros(FT,80),
               F_PPAR_NATURE = zeros(FT,80),
               F_PPAR_BLUE = zeros(FT,80),
               F_PPAR_GREEN = zeros(FT,80),
               F_PPAR_RED = zeros(FT,80));
for dfr in eachrow(df)
    leaf.bio.state.cab = dfr.CHL;
    leaf.bio.state.car = dfr.CHL / 7;
    leaf_spectra!(config, leaf.bio, FT(5));
    dfr.F_APAR_NATURE, dfr.F_PPAR_NATURE = leaf_fapar_fppar(config, leaf, rad_nature);
    dfr.F_APAR_BLUE, dfr.F_PPAR_BLUE = leaf_fapar_fppar(config, leaf, rad_blue);
    dfr.F_APAR_GREEN, dfr.F_PPAR_GREEN = leaf_fapar_fppar(config, leaf, rad_green);
    dfr.F_APAR_RED, dfr.F_PPAR_RED = leaf_fapar_fppar(config, leaf, rad_red);
end;

# 5. save the df to csv
save_csv!("../output/1_absorption.csv", df);
