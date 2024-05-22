#
# this script is used to extend the absorption features of the model to the UV range for the CAB, CAR, and LMA
#
using DataFrames: DataFrame
using LazyArtifacts
using Statistics: mean

using Emerald.EmeraldIO.Text: read_csv, save_csv!
using Emerald.EmeraldMath.Regression: linear_regress
using Emerald.EmeraldLand.Namespace: ReferenceSpectra
using Emerald.EmeraldLand.SPAC: read_spectrum


# 1. read the data
old_spectra = ReferenceSpectra{Float64}(artifact"land_model_spectrum_V5" * "/clima_land_spectra_1nm_2021.nc");
old_df = DataFrame(WL       = old_spectra.Λ,
                   REF_CHL  = old_spectra.K_CAB,
                   REF_CARV = old_spectra.K_CAR_V,
                   REF_CARZ = old_spectra.K_CAR_Z,
                   REF_LMA  = old_spectra.K_LMA,
                   REF_H2O  = old_spectra.K_H₂O);

# 2. extend the spectra for CAB
begin
    # extend the old spectra to match the wavelength bins of the new df
    new_chl = read_csv("../../../data/k_chl.csv");
    new_chl[!,"REF_K"] .= NaN;
    for i in eachindex(new_chl.WL)
        if new_chl.WL[i] >= 400
            new_chl.REF_K[i] = read_spectrum(old_df.WL, old_df.REF_CHL, new_chl.WL[i]);
        end;
    end;

    # rescale the new spectra to match the observed spectra
    y0 = mean(new_chl.REF_K[800 .<= new_chl.WL .<= 1000]);
    y1 = mean(new_chl.REF_K[400 .<= new_chl.WL .<= 402]);
    z0 = mean(new_chl.K_CHL[800 .<= new_chl.WL .<= 1000]);
    z1 = mean(new_chl.K_CHL[400 .<= new_chl.WL .<= 402]);
    slope = (y0 - y1) / (z0 - z1);
    new_chl[!,"PRED_K"] .= NaN;
    for i in eachindex(new_chl.WL)
        if new_chl.WL[i] <= 401
            new_chl.PRED_K[i] = slope * (new_chl.K_CHL[i] - z0) + y0;
        end;
    end;

    # extend the UV part of the chlorophyll absorption spectrum in the old_df file
    for i in eachindex(old_df.WL)
        if old_df.WL[i] < 400
            if old_df.WL[i] > minimum(new_chl.WL)
                old_df.REF_CHL[i] = read_spectrum(new_chl.WL, new_chl.PRED_K, old_df.WL[i]);
            else
                old_df.REF_CHL[i] = new_chl.PRED_K[1];
            end;
        end;
    end;
end;

# 4. per file, extend the UV part of the carotenoid absorption spectrum in the old_df file
files = ["../../../data/k_car-1.csv", "../../../data/k_car-2.csv", "../../../data/k_car-3.csv"];
kcarvs = [];
kcarzs = [];
for file in files
    # extend the old spectra to match the wavelength bins of the new df
    new_car = read_csv(file);
    new_car[!,"REF_KV"] .= NaN;
    new_car[!,"REF_KZ"] .= NaN;
    for i in eachindex(new_car.WL)
        if new_car.WL[i] >= 400
            new_car.REF_KV[i] = read_spectrum(old_df.WL, old_df.REF_CARV, new_car.WL[i]);
            new_car.REF_KZ[i] = read_spectrum(old_df.WL, old_df.REF_CARZ, new_car.WL[i]);
        end;
    end;

    # rescale the new spectra to match the observed spectra
    yv = mean(new_car.REF_KV[400 .<= new_car.WL .<= 402]);
    yz = mean(new_car.REF_KZ[400 .<= new_car.WL .<= 402]);
    z1 = mean(new_car.K_CAR[400 .<= new_car.WL .<= 402]);
    slopev = yv / z1;
    slopez = yz / z1;
    new_car[!,"PRED_KV"] .= NaN;
    new_car[!,"PRED_KZ"] .= NaN;
    for i in eachindex(new_car.WL)
        if new_car.WL[i] <= 401
            new_car.PRED_KV[i] = slopev * new_car.K_CAR[i];
            new_car.PRED_KZ[i] = slopez * new_car.K_CAR[i];
        end;
    end;

    # extend the UV part of the carotenoid absorption spectrum in the old_df file
    kcarv = [];
    kcarz = [];
    for i in eachindex(old_df.WL)
        if old_df.WL[i] < 400
            push!(kcarv, read_spectrum(new_car.WL, new_car.PRED_KV, old_df.WL[i]));
            push!(kcarz, read_spectrum(new_car.WL, new_car.PRED_KZ, old_df.WL[i]));
        end;
    end;
    push!(kcarvs, kcarv);
    push!(kcarzs, kcarz);
end;
for i in eachindex(old_df.WL)
    if old_df.WL[i] < 400
        old_df.REF_CARV[i] = mean([kcarvs[1][i], kcarvs[2][i], kcarvs[3][i]]);
        old_df.REF_CARZ[i] = mean([kcarzs[1][i], kcarzs[2][i], kcarzs[3][i]]);
    end;
end;

# 5. extend the UV part of the LMA absorption spectrum in the old_df file
begin
    # extend the old spectra to match the wavelength bins of the new df
    new_lma = read_csv("../../../data/k_lma.csv");
    new_lma[!,"REF_K"] .= NaN;
    for i in eachindex(new_lma.WL)
        if new_lma.WL[i] >= 400
            new_lma.REF_K[i] = read_spectrum(old_df.WL, old_df.REF_LMA, new_lma.WL[i]);
        end;
    end;

    # rescale the new spectra to match the observed spectra
    # y = 4490.6x - 837.39 R² = 0.998 using the last 5 points
    new_lma[!,"PRED_K"] .= NaN;
    for i in eachindex(new_lma.WL)
        if new_lma.WL[i] <= 401
            new_lma.PRED_K[i] = 4490.6 * new_lma.K_LMA[i] - 837.39;
        end;
    end;

    # extend the UV part of the LMA absorption spectrum in the old_df file
    for i in eachindex(old_df.WL)
        if old_df.WL[i] < 400
            old_df.REF_LMA[i] = read_spectrum(new_lma.WL, new_lma.PRED_K, old_df.WL[i]);
        end;
    end;
end;

# 6. save the rescaled curve
save_csv!(old_df, "../output/1_new_spectra.csv");
