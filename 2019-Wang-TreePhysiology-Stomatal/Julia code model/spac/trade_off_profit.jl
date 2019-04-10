# this function returns the relative profit from trade-off model

function GetTradeoffProfit(parts, envir, up_p, down_p, histories, legacies)
    # result = [soilP, leafP, Ca, E, G, Ci, A, dEdP, vpd, leafTem]
    result = []
    # assign the parameters
    # envir=[tem_a, vpr_a, [CO2], light, speed, r_abs]
    # part[leaf]=[b, c, k, vcmax25, vjratio, leaf_g, leaf/stem, gmax, leaf_w, life_time, mass_per_area]
    tem_a = envir[1]
    vpr_a = envir[2]
    co2_a = envir[3]
    light = envir[4]
    speed = envir[5]
    r_abs = envir[6]
    vcmax25 = parts[4][4]
    vjratio = parts[4][5]
    leaf_g = parts[4][6]
    gmax25 = parts[4][8]
    leaf_w = parts[4][9]
    # get the hydraulic supply cost
    cost = GetTradeoffHydraulicCost(parts, up_p, down_p, histories, legacies)
    emol = cost[3] * leaf_g / gmax25 * 0.015432
    # get leaf temperature and do not reset gmax
    tem_l = GetLeafTemperature(tem_a, speed, leaf_w, emol, r_abs)
    tem_m = 0.5 * (tem_a + tem_l)
    gmax = gmax25
    # get vpd, and do not fix
    vpr_l = GetSaturatedVaporPressure(tem_l)
    vpd = vpr_l - vpr_a
    # get the photosynthetic gain and profit
    if vpd <= 0
        result = [up_p, down_p, co2_a, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    else
        g_tree = cost[3] / vpd * 100.0
        # if g_tree is higher than maximum capability
        if g_tree > gmax
            result = [up_p, down_p, co2_a, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        # whenever possible
        else
            para = [vcmax25, vjratio, tem_l, co2_a]
            gain = GetTradeoffPhotosyntheticGain(para, light, g_tree*leaf_g/gmax25)
            result = [up_p, down_p, co2_a, cost[3], g_tree, gain[2], gain[3], cost[4], vpd, tem_l]
        end
    end
    return result
end
