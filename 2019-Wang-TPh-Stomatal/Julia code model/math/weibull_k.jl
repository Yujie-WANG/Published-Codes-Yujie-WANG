# this function returns the weibull k

function GetWeibullK(b, c, tension)
    if(tension <= 0)
        return 1.0
    else
        return exp( -(tension/b)^c )
    end
end
