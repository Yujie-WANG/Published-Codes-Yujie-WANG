# this function returns the e at given tension

function GetTreeTranspiration(parts, up_p, down_p, histories, legacies)
    # assign soil-plant to parts
    soil = parts[1]
    root = parts[2]
    stem = parts[3]
    leaf = parts[4]
    # define parameters
    min_f = 0.0
    max_f = 1.0
    flow = 0.0
    # if up_p is higher than down_p
    if(down_p <= up_p)
        flow = 0.0
    # if up_p is lower than down_p
    else
        # first find the min and max flow to limit the p
        while true
            tension = GetTreeCanopyPressure(parts, up_p, max_f, histories, legacies)
            if tension > down_p
                break
            else
                min_f = max_f
                max_f *= 2.0
            end
        end
        # use divided method to get flow
        while true
            flow = 0.5 * (max_f + min_f)
            tension = GetTreeCanopyPressure(parts, up_p, flow, histories, legacies)
            if abs(max_f-min_f) < 1E-8
                break
            elseif tension > down_p
                max_f = flow
            else
                min_f = flow
            end
        end
    end
    return flow
end
