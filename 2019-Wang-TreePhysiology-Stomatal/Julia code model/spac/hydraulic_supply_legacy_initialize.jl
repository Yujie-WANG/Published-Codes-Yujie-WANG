# this function update the legacies

function TreeLegaciesInitialize(legacies)
    default = ones(200,2)
    default[:,1] = 0
    legacies[1,:,:] = default
    legacies[2,:,:] = default
    legacies[3,:,:] = default
end
