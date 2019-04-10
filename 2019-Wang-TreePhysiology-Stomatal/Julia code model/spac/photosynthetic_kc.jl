# this function returns the kc from tem

function GetPhotosyntheticKc(tem)
    return CONST_KC * 2.1 ^ (0.1*(tem-25.0))
end
