# This function returns the tissue pressure gradient

function GetPlantDP(tissue, up_p, flow, history, legacy)
    # assign the b c kmax
    b = tissue[1]
    c = tissue[2]
    k = tissue[3]
    # this takes consideration of height impacts
    h = 0.0
    if length(tissue) == 4
        h = tissue[4]
    end
    # operate on each slice
    tension = up_p
    dp = 0.0
    if history == 0
        for i in 1:200
            shell_f = GetWeibullK(b, c, tension)
            shell_k = k * 200.0 * shell_f
            dp += flow / shell_k
            dp += 5E-5 * h
            tension = up_p + dp
        end
    else
        for i = 1:200
            # remember new legacy effect
            if tension > legacy[i][1]
                shell_f = GetWeibullK(b, c, tension)
                legacy[1,i,1] = tension
                legacy[1,i,2] = shell_f
            end
            shell_k = k * 200.0 * legacy[1,i,2]
            dp += flow / shell_k
            dp += 5E-5 * h
            tension = up_p + dp
        end
    end
    return dp
end
