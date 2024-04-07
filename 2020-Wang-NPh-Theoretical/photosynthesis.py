# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 08:57:45 2018

@author: Yujie
"""

from pylab import array,exp,sqrt

# get one-point measurement Vcmax
def GetOnePointVcmax(ci, an, tem, gamma=2.5):
    v_min = 1.0
    v_max = 200.0
    v25   = 100.0
    while v_max-v_min>1:
        v25   = 0.5 * (v_max + v_min)
        r_day = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0)))
        vmax = GetPhotosyntheticVcmax(v25,tem)
        kc = 41.01637 * 2.1**(0.1*(tem-25.0))
        ko = 28201.92 * 1.2**(0.1*(tem-25.0))
        km = kc * (1.0+21000.0/ko)
        ac = vmax * (ci-gamma) / (ci+km)
        af = ac - r_day
        if (af>an and an>0) or (an>af and an<0):
            v_max = v25
        else:
            v_min = v25
    return v25


# calculate j from light
def GetPhotosyntheticJ(jmax, light):
    a = 0.9
    b = -0.3*light - jmax
    c = 0.3*light*jmax
    j = ( -b - sqrt(b*b-4*a*c) ) / a * 0.5
    return j

# calculate jmax from temperature
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

# calculate vcmax from temperature
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

def get_a(v25,j25,gamma,ci,tem,par):
    adjust = 0.98
    tar_p = ci
    r_day  = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0)))
    vmax = GetPhotosyntheticVcmax(v25,tem)
    jmax = GetPhotosyntheticJmax(j25,tem)
    j = GetPhotosyntheticJ(jmax,par)
    kc = 41.01637 * 2.1**(0.1*(tem-25.0))
    ko = 28201.92 * 1.2**(0.1*(tem-25.0))
    km = kc * (1.0+21000.0/ko)
    aj = j * (tar_p-gamma) / (4.0*(tar_p+2*gamma))
    ac = vmax * (tar_p-gamma) / (tar_p+km)
    af = (aj + ac - sqrt((aj+ac)**2.0 - 4*adjust*aj*ac) ) / adjust * 0.5
    af = af - r_day
    return af

def get_a_seg(v25,j25,gamma,ci,tem,par):
    tar_p = ci
    r_day  = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0)))
    vmax = GetPhotosyntheticVcmax(v25,tem)
    jmax = GetPhotosyntheticJmax(j25,tem)
    j = GetPhotosyntheticJ(jmax,par)
    kc = 41.01637 * 2.1**(0.1*(tem-25.0))
    ko = 28201.92 * 1.2**(0.1*(tem-25.0))
    km = kc * (1.0+21000.0/ko)
    aj = j * (tar_p-gamma) / (4.0*(tar_p+2*gamma))
    ac = vmax * (tar_p-gamma) / (tar_p+km)
    af = min(aj,ac)
    af = af - r_day
    return af

# get ci and Anet from gc and ca
def get_a_ci(v25,j25,gamma,gc,ca,tem,par):
    tar_p  = 0.0
    tar_a  = 0.0
    max_p  = ca
    min_p  = gamma
    adjust = 0.98
    r_day  = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0)))
    while(1):
        tar_p = 0.5 * (max_p+min_p)
        vmax = GetPhotosyntheticVcmax(v25,tem)
        jmax = GetPhotosyntheticJmax(j25,tem)
        j = GetPhotosyntheticJ(jmax,par)
        kc = 41.01637 * 2.1**(0.1*(tem-25.0))
        ko = 28201.92 * 1.2**(0.1*(tem-25.0))
        km = kc * (1.0+21000.0/ko)
        aj = j * (tar_p-gamma) / (4.0*(tar_p+2*gamma))
        ac = vmax * (tar_p-gamma) / (tar_p+km)
        af = (aj + ac - sqrt((aj+ac)**2.0 - 4*adjust*aj*ac) ) / adjust * 0.5
        af = af - r_day
        #af = min(aj,ac) - r_day
        tmp_g = af / (ca-tar_p)
        if(abs(tmp_g-gc)/gc < 1E-12):
            tar_a = af
            break
        elif(tmp_g<gc):
            min_p = tar_p
        else:
            max_p = tar_p
        if(abs(max_p-min_p) < 1E-12):
            tar_a = af
            break
    return [tar_p, tar_a]

def get_a_ci_seg(v25,j25,gamma,gc,ca,tem,par):
    tar_p  = 0.0
    tar_a  = 0.0
    max_p  = ca
    min_p  = gamma
    r_day  = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0)))
    while(1):
        tar_p = 0.5 * (max_p+min_p)
        vmax = GetPhotosyntheticVcmax(v25,tem)
        jmax = GetPhotosyntheticJmax(j25,tem)
        j = GetPhotosyntheticJ(jmax,par)
        kc = 41.01637 * 2.1**(0.1*(tem-25.0))
        ko = 28201.92 * 1.2**(0.1*(tem-25.0))
        km = kc * (1.0+21000.0/ko)
        aj = j * (tar_p-gamma) / (4.0*(tar_p+2*gamma))
        ac = vmax * (tar_p-gamma) / (tar_p+km)
        af = min(aj,ac)
        af = af - r_day
        #af = min(aj,ac) - r_day
        tmp_g = af / (ca-tar_p)
        if(abs(tmp_g-gc)/gc < 1E-12):
            tar_a = af
            break
        elif(tmp_g<gc):
            min_p = tar_p
        else:
            max_p = tar_p
        if(abs(max_p-min_p) < 1E-12):
            tar_a = af
            break
    return [tar_p, tar_a]

def get_gamma_nostar(v25,gamma,tem,par):
    tar_p  = 0.0
    max_p  = 80.0
    min_p  = gamma
    r_day  = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0)))
    while(1):
        tar_p = 0.5 * (max_p+min_p)
        vmax = GetPhotosyntheticVcmax(v25,tem)
        kc = 41.01637 * 2.1**(0.1*(tem-25.0))
        ko = 28201.92 * 1.2**(0.1*(tem-25.0))
        km = kc * (1.0+21000.0/ko)
        ac = vmax * (tar_p-gamma) / (tar_p+km)
        af = ac - r_day
        if abs(af)<1E-12:
            break
        elif af>0:
            max_p = tar_p
        else:
            min_p = tar_p
        if(abs(max_p-min_p) < 1E-12):
            break
    return tar_p

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
        rday = v25 * 0.01 * 2**(0.1*(t[i]-25.0)) / (1.0+exp(1.3*(t[i]-55.0)))
        af = max(0, min(aj,ac) ) - rday
        a_list.append(af)
    ym = array(a_list)
    result = ym - y
    return result