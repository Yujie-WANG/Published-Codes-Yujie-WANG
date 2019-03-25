# this function returns leaf temperature form simple version

function GetLeafTemperature(tem_a, wind, leaf_w, emol, rabs)
    # energy from emittance
    e_emit = 0.97 * GetBlackbodyEmittance(tem_a)
    # energy from evaporation
    e_vapo = emol * (45064.3-42.9143*tem_a)
    # boundary condition
    gr = GetRadiativeConductance(tem_a)
    gha = GetBoundaryLayerConductance(wind, leaf_w)
    # leaf temperrature from twosided
    tem_l = tem_a + (rabs-e_emit-0.5*e_vapo) / (29.3*(gr+gha))
    return tem_l
end
