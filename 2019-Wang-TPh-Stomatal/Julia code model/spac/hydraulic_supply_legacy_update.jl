# this function update the legacies

function TreeLegaciesUpdate(parts, up_p, flow, histories, legacies)
    # assign soil-plant to parts
    soil = parts[1]
    root = parts[2]
    stem = parts[3]
    leaf = parts[4]
    tension = up_p
    # assign the legacies
    legacy_r = legacies[1,:,:]
    legacy_s = legacies[2,:,:]
    legacy_l = legacies[3,:,:]
    # get soil dp
    tension += GetSoilDP(soil, tension, flow)
    # get root-stem-leaf dp
    tension += GetPlantDP(root, tension, flow, histories[1], legacy_r)
    tension += GetPlantDP(stem, tension, flow, histories[2], legacy_s)
    tension += GetPlantDP(leaf, tension, flow, histories[3], legacy_l)
    # use history to re-define legacies
    if histories[1] == 0
        legacy_r[:,1] = 0.0
        legacy_r[:,2] = 1.0
    end
    if histories[2] == 0
        legacy_s[:,1] = 0.0
        legacy_s[:,2] = 1.0
    end
    if histories[3] == 0
        legacy_l[:,1] = 0.0
        legacy_l[:,2] = 1.0
    end
    legacies[1,:,:] = legacy_r
    legacies[2,:,:] = legacy_s
    legacies[3,:,:] = legacy_l
end
