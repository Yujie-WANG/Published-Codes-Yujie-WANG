# This function takes input from the STDIN and judge if you have the right input
function GetInputNumber(indent)
    judge = 0
    num_flt = 0.0
    while judge == 0
        if indent == 4
            print_with_color(:yellow, "    Please input a valid number: ")
        elseif indent == 8
            print_with_color(:yellow, "        Please input a valid number: ")
        else
            print_with_color(:yellow, "Please input a valid number: ")
        end
        num_str = strip(readline(STDIN))
        num_flt =
            try
                float(num_str)
            catch
                "The string you input is invalid to convert to float!"
            end
        judge = isa(num_flt, Float64)
    end
    return num_flt
end
