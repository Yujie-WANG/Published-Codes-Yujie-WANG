# this function returns the radiative conductance from simple version

function GetRadiativeConductance(tem)
    cond = 0.1579 + 0.0017*tem + 7.17E-6*tem^2
    return cond
end
