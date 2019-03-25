# this returns the soil water content

function GetSoilRWC(bulkA, bulkN, bulkM, tmpP)
    bulkT = ( 1.0 / (1.0 + (bulkA*tmpP)^bulkN) ) ^ bulkM
    return bulkT
end
