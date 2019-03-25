# this function returns the vcmax from vcmax25 and tem

function GetPhotosyntheticVcmax(vcmax25, tem)
    ha=73637.0
    hd=149252.0
    sv=486.0
    t0=298.15
    r=8.315
    c = 1.0 + exp((sv*t0 -hd)/(r*t0))
    t1 = tem + 273.15
    factor = c * exp(ha/r/t0*(1.0-t0/t1)) / (1.0 + exp((sv*t1-hd)/(r*t1)))
    vcmax = vcmax25 * factor
    return vcmax
end
