# this function returns saturated vapor presure at 0 pressure

function GetSaturatedVaporPressure(tem)
    pres = 0.611 * exp(17.502 * tem / (tem+240.97))
    return pres
end
