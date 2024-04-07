# this function returns the ko from tem

function GetPhotosyntheticKo(tem)
    return CONST_KO * 1.2 ^ (0.1*(tem-25.0))
end
