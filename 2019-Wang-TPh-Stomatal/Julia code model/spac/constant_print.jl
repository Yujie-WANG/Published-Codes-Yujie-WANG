# print the global parameters used in spac model
function PrintSPACConstants()
    print_with_color(:red, "\nCurrent SPAC constants are:\n")
    print_with_color(:black, "    Soil hydraulics --- " * string(PARTS[1][1]) * " " * string(PARTS[1][2]) * " " * string(PARTS[1][3]) * " " * string(PARTS[1][4]) * ";\n")
    print_with_color(:black, "    Root hydraulics --- " * string(PARTS[2][1]) * " " * string(PARTS[2][2]) * " " * string(PARTS[2][3]) * " " * string(PARTS[2][4]) * ";\n")
    print_with_color(:black, "    Stem hydraulics --- " * string(PARTS[3][1]) * " " * string(PARTS[3][2]) * " " * string(PARTS[3][3]) * " " * string(PARTS[3][4]) * ";\n")
    print_with_color(:black, "    Leaf hydraulics --- " * string(PARTS[4][1]) * " " * string(PARTS[4][2]) * " " * string(PARTS[4][3]) * ";\n")
    print_with_color(:black, "    Leaf photosynthesis --- " * string(PARTS[4][4]) * " " * string(PARTS[4][5]) * ";\n")
    print_with_color(:black, "    Leaf gas exchange --- " * string(PARTS[4][6]) * " " * string(PARTS[4][7]) * " " * string(PARTS[4][8]) * ";\n")
    print_with_color(:black, "    Leaf width --- " * string(PARTS[4][9]) * ";\n")
    print_with_color(:black, "    Leaf life time (optimization) --- " * string(PARTS[4][10]) * ";\n")
    print_with_color(:black, "    Leaf mass per area (optimization) --- " * string(PARTS[4][11]) * ";\n")
    print_with_color(:black, "    Leaf hydraulic details (optimization) --- " * string(PARTS[4][12]) * " " * string(PARTS[4][13]) * " " * string(PARTS[4][14]) * ";\n")
end
