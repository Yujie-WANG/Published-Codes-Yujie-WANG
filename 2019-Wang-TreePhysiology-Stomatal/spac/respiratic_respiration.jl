# this function returns respiration rate from vcmax and tem

function GetRespirationRate(vcmax25, tem)
    rday25 = vcmax25 * 0.015
    rday = rday25 * 2.0 ^ (0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0)))
    return rday
end
