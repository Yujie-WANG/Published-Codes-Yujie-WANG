# this function ouput number in given format

function GetOutputNumber(number, digits)
    str_num = repeat(" ", digits)
    
    # transform the number when number > 99999999
    if number > 10^digits-1
        tmp_dig = digits-6
        if tmp_dig == 1
            str_num = @sprintf "%.1e" number
        elseif tmp_dig == 2
            str_num = @sprintf "%.2e" number
        elseif tmp_dig == 3
            str_num = @sprintf "%.3e" number
        elseif tmp_dig == 4
            str_num = @sprintf "%.4e" number
        end
    
    # transform the number when 9999999.5< number < 99999999
    elseif number <= 10^digits-1 && number >= 10^(digits-1)-0.5
        str_num = string( round(Int, number) )
    
    # transform the number when 1000000 <= number < 9999999.5
    elseif number < 10^(digits-1)-0.5 && number >= 10^(digits-2)
        str_num = " " * string( round(Int, number) )
    
    # transform the number when 0 <= number < 1000000
    elseif number >= 0 && number < 10^(digits-2)
        if typeof(number) == Int64
            tmp_str = string(number)
            tmp_spc = repeat(" ", digits-length(tmp_str))
            str_num = tmp_spc * tmp_str
        elseif typeof(number) == Float64
            tmp_str = @sprintf "%.16f" number
            str_num = string( tmp_str[1:digits] )
        end
    
    # transform the number when number < -9999999
    elseif number < -10^(digits-1)+1
        tmp_dig = digits-7
        if tmp_dig == 1
            str_num = @sprintf "%.1e" number
        elseif tmp_dig == 2
            str_num = @sprintf "%.2e" number
        elseif tmp_dig == 3
            str_num = @sprintf "%.3e" number
        elseif tmp_dig == 4
            str_num = @sprintf "%.4e" number
        end
    end
    
    return str_num
end
