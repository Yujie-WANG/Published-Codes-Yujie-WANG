# this returns the soil water pressure

function GetSoilP(bulkA, bulkN, bulkM, bulkT)
    soilP = 1.0 / bulkA * ( bulkT^(-1.0/bulkM) - 1.0 ) ^ (1.0/bulkN)
    return soilP
end
