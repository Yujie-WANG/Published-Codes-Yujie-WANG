from numpy import array,exp,sqrt
from scipy.optimize import leastsq

def GetPhotosyntheticJ(jmax, light):
    a = 0.9
    b = -0.3*light - jmax
    c = 0.3*light*jmax
    j = ( -b - sqrt(b*b-4*a*c) ) / a * 0.5
    return j

def GetPhotosyntheticJmax(jmax25, tem):
    ha=50300.0
    hd=152044.0
    sv=495.0
    t0=298.15
    r=8.315
    c = 1.0 + exp((sv*t0 -hd)/(r*t0))
    t1 = tem + 273.15
    factor = c * exp(ha/r/t0*(1.0-t0/t1)) / (1.0 + exp((sv*t1-hd)/(r*t1)))
    jmax = jmax25 * factor
    return jmax

def GetPhotosyntheticVcmax(vcmax25, tem):
    ha=73637.0
    hd=149252.0
    sv=486.0
    t0=298.15
    r=8.315
    c = 1.0 + exp((sv*t0 -hd)/(r*t0))
    t1 = tem + 273.15
    factor = c * exp(ha/r/t0*(1.0-t0/t1)) / (1.0 + exp((sv*t1-hd)/(r*t1)))
    vcmax = vcmax25 * factor
    return vcmax

def residual_vj(p,x,y,t):
    v25,j25 = p
    a_list = []
    for i in range(len(x)):
        vmax = GetPhotosyntheticVcmax(v25,t[i])
        jmax = GetPhotosyntheticJmax(j25,t[i])
        j = GetPhotosyntheticJ(jmax,1200.0)
        kc = 41.01637 * 2.1**(0.1*(t[i]-25.0))
        ko = 28201.92 * 1.2**(0.1*(t[i]-25.0))
        gamma = 21000.0 * 0.5 / (2600.0*0.57**(0.1*(t[i]-25.0)))
        km = kc * (1.0+21000.0/ko)
        aj = j * (x[i]-gamma) / (4.0*(x[i]+2*gamma))
        ac = vmax * (x[i]-gamma) / (x[i]+km)
        if(x[i] < gamma):
            aj = 0.0
            ac = 0.0
        rday = 1.95 * 2**(0.1*(t[i]-25.0)) / (1.0+exp(1.3*(t[i]-55.0)))
        af = min(aj,ac)# - rday
        a_list.append(af)
    ym = array(a_list)
    result = ym - y
    for i in range(len(x)):
        print(x[i], "\t", y[i], "\t", ym[i], "\t", t[i])
    print("")
    return result

def residual_vjgamma(p,x,y,t):
    v25,j25,gamma = p
    a_list = []
    for i in range(len(x)):
        vmax = GetPhotosyntheticVcmax(v25,t[i])
        jmax = GetPhotosyntheticJmax(j25,t[i])
        j = GetPhotosyntheticJ(jmax,1200.0)
        kc = 41.01637 * 2.1**(0.1*(t[i]-25.0))
        ko = 28201.92 * 1.2**(0.1*(t[i]-25.0))
        #gamma = 21000.0 * 0.5 / (2600.0*0.57**(0.1*(t[i]-25.0)))
        km = kc * (1.0+21000.0/ko)
        aj = j * (x[i]-gamma) / (4.0*(x[i]+2*gamma))
        ac = vmax * (x[i]-gamma) / (x[i]+km)
        # 1.95 is the rday for water birch, change it accordingly
        rday = 1.95 * 2**(0.1*(t[i]-25.0)) / (1.0+exp(1.3*(t[i]-55.0)))
        af = max(0, min(aj,ac) ) - rday
        #af = min(aj,ac)
        a_list.append(af)
    ym = array(a_list)
    result = ym - y
    error = sum(result**2)
    for i in range(len(x)):
        print(x[i], "\t", y[i], "\t", ym[i], "\t", t[i])
    print("\n", v25, "\t", j25, "\t", gamma, "\t", error)
    print("")
    return result

def aci_curve(x = range(10), y = range(10), t = range(10), rday=True):
    xx = array(x)/10.0
    yy = array(y)
    tt = array(t)
    if(rday == True):
        p0 = [50.0,100.0]
        plsq = leastsq(residual_vj, p0, args=(xx,yy,tt))
        return [plsq[0][0], plsq[0][1]]
    else:
        p0 = [50.0,80.0,6.0]
        plsq = leastsq(residual_vjgamma, p0, args=(xx,yy,tt))
        return [plsq[0][0], plsq[0][1], plsq[0][2]]
