# this function returns the photosynthesis rate

function GetPhotosyntheticRate(para, light, ci)
    # assign the parameters
    vcmax25 = para[1]
    ratio = para[2]
    tem = para[3]
    jmax25 = vcmax25 * ratio
    vcmax = GetPhotosyntheticVcmax(vcmax25, tem)
    jmax = GetPhotosyntheticJmax(jmax25, tem)
    # 1. get j from jmax and light
    j = GetPhotosyntheticJ(jmax, light)
    # 2. get co2 limited process
    kc = GetPhotosyntheticKc(tem)
    ko = GetPhotosyntheticKo(tem)
    gamma = GetPhotosyntheticGamma(tem)
    km = kc * (1.0 + CONST_PO/ko)
    # 3. get aj, ac, and af
    adjust = 0.98
    aj = j * (ci-gamma) / (4.0*(ci+2*gamma))
    ac = vcmax * (ci-gamma) / (ci+km)
    af = (aj + ac - sqrt((aj+ac)^2 - 4*adjust*aj*ac) ) / adjust * 0.5
    return af
end
