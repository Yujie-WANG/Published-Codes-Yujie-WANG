# this function returns jmax from jmax25 and temperature

function GetPhotosyntheticJmax(jmax25, tem)
    ha=50300.0
    hd=152044.0
    sv=495.0
    t0=298.15
    r=8.315
    c = 1.0 + exp((sv*t0 -hd)/(r*t0))
    t1 = tem + 273.15
    factor = c * exp(ha/r/t0*(1.0-t0/t1)) / (1.0 + exp((sv*t1-hd)/(r*t1)))
    jmax = jmax25 * factor
    return jmax
end
