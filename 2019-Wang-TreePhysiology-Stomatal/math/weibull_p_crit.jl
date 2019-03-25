# this function returns the p_crit

function GetWeibullPCrit(b, c)
    return b * log(1000.0)^(1/c)
end
