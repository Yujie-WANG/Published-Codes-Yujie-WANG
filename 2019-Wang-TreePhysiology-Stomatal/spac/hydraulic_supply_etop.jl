# this function returns the distal tension of the plant

function GetTreeCanopyPressure(parts, up_p, flow, histories, legacies)
    # assign soil-plant to parts
    soil = parts[1]
    root = parts[2]
    stem = parts[3]
    leaf = parts[4]
    tension = up_p
    # get soil dp
    tension += GetSoilDP(soil, tension, flow)
    # get root-stem-leaf dp
    # note that here legacies will not be modified
    tension += GetPlantDP(root, tension, flow, histories[1], legacies[1,:,:])
    tension += GetPlantDP(stem, tension, flow, histories[2], legacies[2,:,:])
    tension += GetPlantDP(leaf, tension, flow, histories[3], legacies[3,:,:])
    return tension
end
