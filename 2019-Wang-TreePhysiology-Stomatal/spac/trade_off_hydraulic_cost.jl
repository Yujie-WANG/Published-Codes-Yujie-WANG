# this function returns the hydraulic cost in trade-off model

function GetTradeoffHydraulicCost(parts, up_p, down_p, histories, legacies)
    # please note that histories will not be saved here
    dp = 1E-3
    f = GetTreeTranspiration(parts, up_p, down_p, histories, legacies)
    f_dp = GetTreeTranspiration(parts, up_p, down_p+dp, histories, legacies)
    s = (f_dp-f) / dp
    # result format in [P, P, E, dEdP]
    return [up_p, down_p, f, s]
end
