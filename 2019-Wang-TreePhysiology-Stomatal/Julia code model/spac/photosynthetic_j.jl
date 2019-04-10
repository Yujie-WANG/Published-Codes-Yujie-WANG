# this function returns the j from jmax and light

function GetPhotosyntheticJ(jmax, light)
    a = CONST_J_THETA
    b = -CONST_J_ALPHA*light - jmax
    c = CONST_J_ALPHA * light * jmax
    j = ( -b - sqrt(b^2-4*a*c) ) / a * 0.5
    return j
end
