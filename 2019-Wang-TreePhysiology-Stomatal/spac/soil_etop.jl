# This function returns the soil pressure gradient

function GetSoilDP(soil, up_p, flow)
    # assign the soil parameters
    a = soil[1]
    n = soil[2]
    m = soil[3]
    k = soil[4]
    tension = up_p
    # iterate through each shell to get dp
    dp = 0.0
    for i = 0:99
        shell_t = (1.0 / (1.0 + (a*tension)^n)) ^ m;
        shell_f = sqrt(shell_t) * (1.0 - (1.0-shell_t^(1.0/m)) ^ m) ^ 2.0
        shell_k = k * shell_f * log(10.0) / log((10.0-0.09*i)/(10.0-0.09*(i+1)))
        dp += flow / shell_k
        tension = up_p + dp
    end
    return dp
end
