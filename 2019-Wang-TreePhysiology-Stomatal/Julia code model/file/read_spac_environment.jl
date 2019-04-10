# this function is meant to read the environmental factors of spac fitting
function ReadSPACEnvironment(file_name)
    file_read = open(file_name, "r")
    while !eof(file_read)
        tmp = readline(file_read)
        if tmp[1]>='0' && tmp[1] <='9'
            list = float( split(tmp,"\t") )
            push!(ARRAY_TAIR , list[1])
            push!(ARRAY_TLEAF, list[2])
            push!(ARRAY_VPA  , list[3])
            push!(ARRAY_CO2  , list[4])
            push!(ARRAY_PAR  , list[5])
            push!(ARRAY_WIND , list[6])
            push!(ARRAY_RABS , list[7])
            push!(ARRAY_PSOIL, list[8])
            push!(ARRAY_AREA , list[9])
            push!(ARRAY_LABA , list[10])
            push!(ARRAY_VCMAX, list[11])
            push!(ARRAY_JMAX , list[12])
            push!(ARRAY_PHOTO, list[14])
            push!(ARRAY_TRANS, list[15])
            push!(ARRAY_PLEAF, list[16])
            push!(ARRAY_KTREE, list[17])
        end
    end
    close(file_read)
end
