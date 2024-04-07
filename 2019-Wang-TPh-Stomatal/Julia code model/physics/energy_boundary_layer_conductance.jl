# this function returns the boundary layer conductance

function GetBoundaryLayerConductance(wind, leaf_w)
    cond = 1.4 * 0.135 * sqrt( wind/(0.72*leaf_w) )
    return cond
end
