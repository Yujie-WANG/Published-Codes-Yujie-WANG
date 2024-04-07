# this function returns the emittance of black body

function GetBlackbodyEmittance(tem)
    kelvin = tem + 273.15
    emit = BOLTZMANN_CONSTANT * kelvin^4
    return emit
end
