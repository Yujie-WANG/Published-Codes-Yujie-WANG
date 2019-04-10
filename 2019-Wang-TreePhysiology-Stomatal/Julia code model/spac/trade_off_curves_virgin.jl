# This function generate the curves so as to draw the virgin curves
function MakeTradeoffCurvesVirgin(up_p)
    virgin_history = [0,0,0]
    virgin_legacies =  ones(3,200,2)
    virgin_legacies[:,:,1] = 0.0
    print_with_color(:red, "\nThe virgin trade-off curves at -" * string(up_p) * " MPa is listed below:\n")
    println("CanopyTension\tFlow\tPhotosynthesis\tdE/DP")
    # iterate through the p range
    down_p = up_p
    crit_p = GetWeibullPCrit(PARTS[4][1], PARTS[4][2])
    dp = 0.01 * (crit_p-up_p)
    while(down_p <=6.0)
        result = GetTradeoffProfitFix(PARTS, ENVIR, up_p, down_p, virgin_history, virgin_legacies, ESFACTOR)
        println(result[2], " ", result[4], " ", result[7], " ", result[8])
        down_p += dp
    end
    print_with_color(:red, "Finished the virgin trade-off curves!\n")
end
