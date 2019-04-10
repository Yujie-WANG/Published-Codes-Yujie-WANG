# This function edit soil type with most known soil types
function EditSoilType()
    # list available soil types
    print_with_color(:red, "\nPlease select the soil types below:\n")
    print_with_color(:blue, "    01. Sand              --- (1479.5945, 2.68)\n")
    print_with_color(:blue, "    02. Loamy Sand        --- (1265.3084, 2.28)\n")
    print_with_color(:blue, "    03. Sandy Loam        --- ( 765.3075, 1.89)\n")
    print_with_color(:blue, "    04. Loamy Sand 2      --- ( 367.3476, 1.56)\n")
    print_with_color(:blue, "    05. Silt              --- ( 163.2656, 1.37)\n")
    print_with_color(:blue, "    06. Silty Loam        --- ( 204.0820, 1.41)\n")
    print_with_color(:blue, "    07. Sandy Clay Loam   --- ( 602.0419, 1.48)\n")
    print_with_color(:blue, "    08. Clay Loam         --- ( 193.8779, 1.31)\n")
    print_with_color(:blue, "    09. Silt Clay Loam    --- ( 102.0410, 1.23)\n")
    print_with_color(:blue, "    10. Sandy Clay Loam 2 --- ( 275.5107, 1.23)\n")
    print_with_color(:blue, "    11. Silty Clay        --- (  51.0205, 1.09)\n")
    print_with_color(:blue, "    12. Clay              --- (  81.6328, 1.09)\n")
    print_with_color(:blue, "    Default. Return to previous menu\n")
    # make the choice
    print_with_color(:red, "Please indicate your choice: ")
    choice = strip(readline(STDIN))
    if choice == "01"
        PARTS[1][1] = 1479.5945;
        PARTS[1][2] = 2.68;
    elseif choice == "02"
        PARTS[1][1] = 1265.3084;
        PARTS[1][2] = 2.28;
    elseif choice == "03"
        PARTS[1][1] = 765.3075;
        PARTS[1][2] = 1.89;
    elseif choice == "04"
        PARTS[1][1] = 367.3476;
        PARTS[1][2] = 1.56;
    elseif choice == "05"
        PARTS[1][1] = 163.2656;
        PARTS[1][2] = 1.37;
    elseif choice == "06"
        PARTS[1][1] = 204.0820;
        PARTS[1][2] = 1.41;
    elseif choice == "07"
        PARTS[1][1] = 602.0419;
        PARTS[1][2] = 1.48;
    elseif choice == "08"
        PARTS[1][1] = 193.8779;
        PARTS[1][2] = 1.31;
    elseif choice == "09"
        PARTS[1][1] = 102.0410;
        PARTS[1][2] = 1.23;
    elseif choice == "10"
        PARTS[1][1] = 275.5107;
        PARTS[1][2] = 1.23;
    elseif choice == "11"
        PARTS[1][1] = 51.0205;
        PARTS[1][2] = 1.09;
    elseif choice == "12"
        PARTS[1][1] = 81.6328;
        PARTS[1][2] = 1.09;
    else
        print_with_color(:red, "    No changes applied!\n")
    end
    PARTS[1][3] = 1.0 - 1.0/PARTS[1][2];
    print_with_color(:red, "    Please check the details to make sure you selected the right soil type!\n")
end
