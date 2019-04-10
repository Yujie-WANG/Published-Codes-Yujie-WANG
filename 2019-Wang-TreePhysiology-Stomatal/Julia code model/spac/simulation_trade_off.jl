# this function is meant to run simulations while you know the environemtal parameters

function SimulateTradeoff(histories, legacies, esfactor)
    # initialize the arrays
    global ARRAY_TAIR  = []
    global ARRAY_TLEAF = []
    global ARRAY_VPA   = []
    global ARRAY_CO2   = []
    global ARRAY_PAR   = []
    global ARRAY_WIND  = []
    global ARRAY_RABS  = []
    global ARRAY_PSOIL = []
    global ARRAY_AREA  = []
    global ARRAY_LABA  = []
    global ARRAY_VCMAX = []
    global ARRAY_JMAX  = []
    global ARRAY_PHOTO = []
    global ARRAY_TRANS = []
    global ARRAY_PLEAF = []
    global ARRAY_KTREE = []

    # assign the value to bo2017
    global parts = ConstantLoadBO2017()

    # read the data sets
    ReadSPACEnvironment("./data/data_file.txt")
    println("Reading finished!")

    # remember the array_psoil and modify it
    MEM_PSOIL = deepcopy(ARRAY_PSOIL)
    FixSimulationKsoil(histories, legacies, esfactor)
    tmp_error = GetSimulationErrorKnowTleaf(histories, legacies, esfactor)

end

function FixSimulationKsoil(histories, legacies, esfactor)
    # run the simulations
    klist = (1:30) * 1E12
    elist = []
    for k in klist
        parts[1][4] = k
        tmp_error = GetSimulationErrorKnowTleaf(histories, legacies, esfactor)
        #print(legacies, "\n")
        push!(elist, tmp_error)
        if length(elist)>=3
            if minimum(elist)==elist[length(elist)-1] && elist[length(elist)-1]<elist[length(elist)] && elist[length(elist)-1]<elist[length(elist)-2]
                break
            end
        end
    end

    # set up best ksoil
    parts[1][4] = klist[length(elist)-1]
    for i in 1:length(elist)
        print(klist[i], "\t", elist[i], "\n")
    end
end

function GetSimulationErrorKnowTleaf(histories, legacies, esfactor)
    # run the simulations
    array_tleaf = []
    array_photo = []
    array_trans = []
    array_pleaf = []
    print("\nSoilP\tLeafP\tE\tG\tCi\tA\tVPD\tLeafT\n")
    for i in 1:length(ARRAY_TAIR)
        # load the environmental factor into model
        parts[2][3] = ARRAY_KTREE[i] * 1.863013458
        parts[3][3] = ARRAY_KTREE[i] * 4.11942096
        parts[4][3] = ARRAY_KTREE[i] * 4.535504046
        parts[4][4] = ARRAY_VCMAX[i]
        parts[4][5] = ARRAY_JMAX[i] / ARRAY_VCMAX[i]
        parts[4][7] = ARRAY_LABA[i]
        parts[4][8] = parts[4][6] * parts[4][7]
        global envir = [ARRAY_TAIR[i], ARRAY_VPA[i], ARRAY_CO2[i], ARRAY_PAR[i], ARRAY_WIND[i], ARRAY_RABS[i]]
        global up_p = ARRAY_PSOIL[i]
        # get the result of trade-off
        # result = [soilP, leafP, Ca, E, G, Ci, A, dEdP, vpd, leafTem]
        result = GetTradeoffProfitOptimizeFixKnowTleaf(parts, envir, up_p, histories, legacies, esfactor, 0, ARRAY_TLEAF[i])
        push!(array_tleaf, result[10])
        push!(array_photo, result[7])
        push!(array_trans, result[4]*ARRAY_AREA[i]/3.6)
        push!(array_pleaf, result[2])
        print(ARRAY_PSOIL[i], "\t", result[2], "\t", result[4], "\t", result[5], "\t", result[6], "\t", result[7], "\t", result[9], "\t", result[10], "\n")
    end

    # calculate the normalized squaure sum of four fittings
    #error_tleaf = array_tleaf - ARRAY_TLEAF
    error_photo = array_photo - ARRAY_PHOTO
    error_trans = array_trans - ARRAY_TRANS
    error_pleaf = array_pleaf - ARRAY_PLEAF
    for i in 1:length(error_pleaf)
        if ARRAY_PLEAF[i] == -1.0
            error_pleaf[i] = NaN
        end
    end
    #std_tleaf = std( convert(Array{Float64,1},ARRAY_TLEAF) )
    std_photo = std( convert(Array{Float64,1},ARRAY_PHOTO) )
    std_trans = std( convert(Array{Float64,1},ARRAY_TRANS) )
    std_pleaf = std( convert(Array{Float64,1},ARRAY_PLEAF[!isnan( convert(Array{Float64,1},error_pleaf) ) ]) )
    sum_error = 0.0
    #sum_error += sum( (error_tleaf/std_tleaf) .^ 2.0 )
    sum_error += sum( (error_photo/std_photo) .^ 2.0 )
    sum_error += sum( (error_trans/std_trans) .^ 2.0 )
    sum_error += sum( (error_pleaf[!isnan( convert(Array{Float64,1},error_pleaf) )]/std_pleaf) .^ 2.0 )
    print(sum_error, "\n")
    return sum_error
end

