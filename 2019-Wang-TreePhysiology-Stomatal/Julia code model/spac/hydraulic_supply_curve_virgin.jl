# This function provides the supply curve at current soil water potential
function MakeHydraulicSupplyCurveVirgin(up_p)
    virgin_history = [0,0,0]
    virgin_legacies =  ones(3,200,2)
    virgin_legacies[:,:,1] = 0.0
    print_with_color(:red, "\nThe virgin supply curve at current soil moisture is listed below:\n")
    println("CanopyTension\tFlow")
    # iterate through the p range
    down_p = up_p
    crit_p = GetWeibullPCrit(PARTS[4][1], PARTS[4][2])
    dp = 0.01 * (crit_p-up_p)
    while(down_p <= crit_p)
        flow = GetTreeTranspiration(PARTS, up_p, down_p, virgin_history, virgin_legacies)
        println(down_p, "\t", flow)
        down_p += dp
    end
    print_with_color(:red, "Finished the virgin supply curve!\n")
end
