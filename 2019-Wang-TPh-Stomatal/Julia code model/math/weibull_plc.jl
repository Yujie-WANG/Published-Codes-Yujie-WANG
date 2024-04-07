# this function returns the weibull plc

function GetWeibullPLC(b, c, tension)
    k = GetWeibullK(b, c, tension)
    return 100 * (1-k)
end
