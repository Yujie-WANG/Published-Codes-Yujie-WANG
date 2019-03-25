# edit the global parameters defined in constant.jl
function EditSPACConstants()
    judge = 1
    while judge == 1
        # make the choice
        choice = MakeSPACConstantsChoice()
        if choice == "00"
            PrintSPACConstants()
        elseif choice == "01"
            EditSoilType()
            print_with_color(:red, "    Please input soil k below\n")
            PARTS[1][4] = GetInputNumber(8)
        elseif choice == "02"
            print_with_color(:red, "    Please input root b below\n")
            PARTS[2][1] = GetInputNumber(8)
            print_with_color(:red, "    Please input root c below\n")
            PARTS[2][2] = GetInputNumber(8)
            print_with_color(:red, "    Please input root k below\n")
            PARTS[2][3] = GetInputNumber(8)
            print_with_color(:red, "    Please input root depth (used for gravity correction) below\n")
            PARTS[2][4] = GetInputNumber(8)
        elseif choice == "03"
            print_with_color(:red, "    Please input stem b below\n")
            PARTS[3][1] = GetInputNumber(8)
            print_with_color(:red, "    Please input stem c below\n")
            PARTS[3][2] = GetInputNumber(8)
            print_with_color(:red, "    Please input stem k below\n")
            PARTS[3][3] = GetInputNumber(8)
            print_with_color(:red, "    Please input stem height (used for gravity correction) below\n")
            PARTS[3][4] = GetInputNumber(8)
        elseif choice == "04"
            print_with_color(:red, "    Please input leaf b below\n")
            PARTS[4][1] = GetInputNumber(8)
            print_with_color(:red, "    Please input leaf c below\n")
            PARTS[4][2] = GetInputNumber(8)
            print_with_color(:red, "    Please input leaf k below\n")
            PARTS[4][3] = GetInputNumber(8)
        elseif choice == "05"
            print_with_color(:red, "    Please input leaf Vcmax25 below\n")
            PARTS[4][4] = GetInputNumber(8)
            print_with_color(:red, "    Please input leaf Vcmax25:Jmax25 ratio below\n")
            PARTS[4][5] = GetInputNumber(8)
        elseif choice == "06"
            print_with_color(:red, "    Please input G per leaf area below\n")
            PARTS[4][6] = GetInputNumber(8)
            print_with_color(:red, "    Please input leaf:stem ratio below\n")
            PARTS[4][7] = GetInputNumber(8)
            print_with_color(:red, "    Please input Gmax per basal area below\n")
            PARTS[4][8] = GetInputNumber(8)
        elseif choice == "07"
            print_with_color(:red, "    Please input leaf width (used for boundary layer) below\n")
            PARTS[4][9] = GetInputNumber(8)
        elseif choice == "08"
            print_with_color(:red, "    Please input leaf life time below\n")
            PARTS[4][10] = GetInputNumber(8)
        elseif choice == "09"
            print_with_color(:red, "    Please input leaf mass per leaf area below\n")
            PARTS[4][11] = GetInputNumber(8)
        elseif choice == "10"
            print_with_color(:yellow, "The function has yet to be added!\n")
        else
            judge = 0
        end
    end
    print_with_color(:red, "\nFinish editting the constants\n")
end

function MakeSPACConstantsChoice()
    # print the choices
    print_with_color(:red, "\nPlease indicate which parameter(s) you want to change:\n")
    print_with_color(:blue, "    00. Print current settings;\n")
    print_with_color(:blue, "    01. Soil hydraulics (alpha, n, m, ksoil);\n")
    print_with_color(:blue, "    02. Root hydraulics (b, c, kmax, depth);\n")
    print_with_color(:blue, "    03. Stem hydraulics (b, c, kmax, height);\n")
    print_with_color(:blue, "    04. Leaf hydraulics (b, c, kmax);\n")
    print_with_color(:blue, "    05. Leaf photosynthesis (vcmax25, vjratio);\n")
    print_with_color(:blue, "    06. Leaf gas exchange (leaf_g, ls, gmax);\n")
    print_with_color(:blue, "    07. Leaf width (leaf_w);\n")
    print_with_color(:blue, "    08. Optimization --- Leaf life time (leaf_t);\n")
    print_with_color(:blue, "    09. Optimization --- Leaf mass per area (lma);\n")
    print_with_color(:blue, "    10. Optimization --- Leaf hydraulic details (leaf_s, leaf_n, leaf_k);\n")
    print_with_color(:blue, "    Default. Back to main menu\n")
    # select which to change
    print_with_color(:red, "Please indicate your choice: ")
    choice = strip(readline(STDIN))
    return choice
end
