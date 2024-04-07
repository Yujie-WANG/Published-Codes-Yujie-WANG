# this function optimize profit by altering leaf water potential

function GetTradeoffProfitOptimizeFixKnowTleaf(parts, envir, up_p, histories, legacies, esfactor, print, tleaf)
    # generate first three components of the chain
    curve = []
    crit_p = GetWeibullPCrit(parts[4][1], parts[4][2])
    down_p = up_p
    result = GetTradeoffProfitFixKnowTleaf(parts, envir, up_p, down_p, histories, legacies, esfactor, tleaf)
    push!(curve, result)
    down_p = crit_p
    result = GetTradeoffProfitFixKnowTleaf(parts, envir, up_p, down_p, histories, legacies, esfactor, tleaf)
    push!(curve, result)
    down_p = 0.5 * (crit_p+up_p)
    result = GetTradeoffProfitFixKnowTleaf(parts, envir, up_p, down_p, histories, legacies, esfactor, tleaf)
    push!(curve, result)
    find(x->( x[4]==0 ), curve)
    # use di-section method to find physical limitation
    while true
        # sort by down_p
        sort!(curve, by=x->(x[2]))
        tmp = find(x->( x[4]==0 ), curve)
        if length(tmp) > 1
            next_site = tmp[2]
            if curve[next_site][2]-curve[next_site-1][2] < 1E-6
                break
            end
            down_p = (curve[next_site][2] + curve[next_site-1][2]) * 0.5
            result = GetTradeoffProfitFixKnowTleaf(parts, envir, up_p, down_p, histories, legacies, esfactor, tleaf)
            push!(curve, result)
        else
            break
        end
    end
    # use di-section method to find optimal
    opt_site = 0
    while true
        # 0. sort the curve by P
        sort!(curve, by=x->(x[2]))
        # 1. sort the tmp_curve by A find the maximal A
        tmp_curve = sort(curve, by=x->(x[7]), rev=true)
        opt_a = tmp_curve[1][7]
        # 2. compute the profit along the list
        opt_site = 1
        opt_prof = 0
        for i in 2:(length(curve)-1)
            cost_r = 1.0 - curve[i][8]/curve[1][8]
            gain_r = curve[i][7] / opt_a
            prof_r = gain_r - cost_r
            if prof_r > opt_prof
                opt_prof = prof_r
                opt_site = i
            end
        end
        # 3. make judgement if break
        if opt_site > 1 && opt_site < length(curve)
            if curve[opt_site+1][2] - curve[opt_site-1][2] < 1E-6
                break
            end
        elseif opt_site == length(curve)
            if curve[opt_site][2] - curve[opt_site-1][2] < 1E-6
                break
            end
        else
            if curve[opt_site+1][2] - curve[opt_site][2] < 1E-6
                break
            end
        end
        # 4. add left point to curve
        if opt_site > 1
            down_p = (curve[opt_site][2] + curve[opt_site-1][2]) * 0.5
            result = GetTradeoffProfitFixKnowTleaf(parts, envir, up_p, down_p, histories, legacies, esfactor, tleaf)
            push!(curve, result)
        end
        # 5. add right point to curve
        if opt_site < length(curve)
            down_p = (curve[opt_site][2] + curve[opt_site+1][2]) * 0.5
            result = GetTradeoffProfitFixKnowTleaf(parts, envir, up_p, down_p, histories, legacies, esfactor, tleaf)
            push!(curve, result)
        end
    end
    # print ouy curves
    if(print == "y")
        for i in curve
            println(i)
        end
    end
    # now opt site is know, result can be returned
    return curve[opt_site]
end
