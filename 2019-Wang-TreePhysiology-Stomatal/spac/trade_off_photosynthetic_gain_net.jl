# this function returns the photosynthetic gain at given conditions

function GetTradeoffPhotosyntheticGainNet(para, light, gw)
    # here para is [vcmax25, vjratio, tem, [CO2]]
    gamma = GetPhotosyntheticGamma(para[3])
    pca = para[4]
    # gw here is in Kg m-2 h-1, 1000/18/3600*1E-5(into per Pa)*1E6(into umol)
    tar_g = gw * 0.154321 / 1.6
    tar_p = 0.0
    tar_a = 0.0
    # set max_p, min_p and tar_p
    max_p = pca
    min_p = gamma
    while true
        tar_p = 0.5 * (max_p+min_p)
        tar_r = GetRespirationRate(para[1], para[3])
        tar_a = GetPhotosyntheticRate(para, light, tar_p) - tar_r
        tmp_g = tar_a / (pca-tar_p)
        if abs(tmp_g-tar_g) < 1E-10
            break
        elseif tmp_g < tar_g
            min_p = tar_p
        else
            max_p = tar_p
        end
    end
    return [pca, tar_p, tar_a]
end
