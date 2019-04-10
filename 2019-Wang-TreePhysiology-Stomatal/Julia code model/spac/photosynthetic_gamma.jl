# this function returns tao from tem

function GetPhotosyntheticGamma(tem)
    t = 2600.0 * 0.57 ^ (0.1*(tem-25.0))
    # gamma = 0.5 * CONST_PO / t
    gamma = 1.8
    return gamma
end
