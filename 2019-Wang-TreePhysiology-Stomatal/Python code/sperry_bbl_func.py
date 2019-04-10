# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 08:49:42 2018

@author: Yujie
"""

# import path and add new path to the list
from sys import path,platform
if platform=="win32":
    cust_folder = "C:\\Users\\Yujie\\Documents\\MEGA\\Program\\yujie_class\\"
else:
    cust_folder = "/home/yujie/MEGA/Program/yujie_class/"
if cust_folder not in path:
    path.append(cust_folder)

# import modules
from numpy import array,exp,log,nanmean,nanstd,sqrt

from photosynthesis  import get_a

# define the function to get p_leaf
def get_p_leaf_birch(ps, km, e):
    rhiz_k = km * 6E9
    rhiz_a = 367.3476
    rhiz_n = 1.56
    rhiz_m = 1.0 - 1.0/rhiz_n
    root_b = 1.879
    root_c = 2.396
    root_k = km * 1.863013458
    stem_b = 2.238
    stem_c = 9.380
    stem_k = km * 4.11942096
    leaf_b = 1.897
    leaf_c = 2.203
    leaf_k = km * 4.535504046
    up_p   = ps
    down_p = ps
    dp     = 0.0
    for i in range(10):
        if down_p>0:
            shell_t = (1.0 / (1.0 + (rhiz_a*down_p)**rhiz_n)) ** rhiz_m
        else:
            shell_t = 1.0
        shell_f  = sqrt(shell_t) * (1.0 - (1.0-shell_t**(1.0/rhiz_m)) ** rhiz_m) ** 2.0
        shell_k  = max( rhiz_k * shell_f * log(10.0) / log((10.0-0.9*i)/(10.0-0.9*(i+1))), 1E-12 )
        dp      += e/shell_k
        down_p   = up_p + dp
    for i in range(100):
        tmp_k = 100.0 * root_k * exp(-(down_p/root_b)**root_c)
        down_p += e / tmp_k
    for i in range(100):
        tmp_k = 100.0 * stem_k * exp(-(down_p/stem_b)**stem_c)
        down_p += e / tmp_k
    for i in range(100):
        tmp_k = 100.0 * leaf_k * exp(-(down_p/leaf_b)**leaf_c)
        down_p += e / tmp_k
    return down_p

# define the function to obtain cross point of empirical function and A-Ci
def get_bbl_prediction(p, vm, jm, ca, d, ps, tem):
    ci_max = ca
    ci_min = 2.5
    while True:
        ci = 0.5 * (ci_max+ci_min)
        a  = get_a(vm,jm,2.5,ci,tem,1000.0)
        gc = a / (ca-ci) * 1E-6 * 101315.0
        gh = 1.6 * gc
        gs = p[0] + 1.6*a/ca * 0.1 * (1.0+p[1]/sqrt(0.01*d)) * (p[2]-ps)/(p[2]-p[3])
        if abs(gs-gh)<1E-5 or abs(ci_max-ci_min)<1E-5:
            break
        elif gs>gh:
            ci_min = ci
        else:
            ci_max = ci
    return a,gs

# define the residual function
def residual_bbl(p, data_tmp):
    print(p)
    print("\n")
    list_d   = data_tmp["D"   ].values
    list_a   = data_tmp["A"   ].values
    list_ca  = data_tmp["Ca"  ].values
    list_ps  = data_tmp["Ps"  ].values
    list_vm  = data_tmp["V"   ].values
    list_jm  = data_tmp["J"   ].values
    list_lb  = data_tmp["LaBa"].values
    list_pl  = data_tmp["Pl"  ].values
    list_e   = data_tmp["E"   ].values
    list_t   = data_tmp["T"   ].values
    list_km  = data_tmp["Km"  ].values
    subdata  = data_tmp.query("Pl>0")
    std_a    = nanstd(list_a       )
    std_e    = nanstd(list_e       )
    std_p    = nanstd(subdata["Pl"])
    residual = []
    for i in range(len(list_d)):
        a,g = get_bbl_prediction(p, list_vm[i], list_jm[i], list_ca[i], list_d[i], list_ps[i], list_t[i])
        e   = g * list_d[i] * 0.648 * list_lb[i] # 0.648 convert this to kg/(h m2)
        pl  = get_p_leaf_birch(list_ps[i], list_km[i], e)
        residual.append( sqrt( abs(list_a[i]-a)/std_a ) )
        residual.append( sqrt( abs(list_e[i]-e)/std_e ) )
        if list_pl[i]>0:
            residual.append( sqrt( abs(list_pl[i]-pl)/std_p ) )
    return array(residual)

def pred_bbl(p, data_tmp):
    list_d   = data_tmp["D"   ].values
    list_ca  = data_tmp["Ca"  ].values
    list_ps  = data_tmp["Ps"  ].values
    list_vm  = data_tmp["V"   ].values
    list_jm  = data_tmp["J"   ].values
    list_lb  = data_tmp["LaBa"].values
    list_t   = data_tmp["T"   ].values
    list_km  = data_tmp["Km"  ].values
    pred_a   = []
    pred_e   = []
    pred_p   = []
    pred_g   = []
    for i in range(len(list_d)):
        a,g = get_bbl_prediction(p, list_vm[i], list_jm[i], list_ca[i], list_d[i], list_ps[i], list_t[i])
        e   = g * list_d[i] * 0.648 * list_lb[i] # 0.648 convert this to kg/(h m2)
        pl  = get_p_leaf_birch(list_ps[i], list_km[i], e)
        pred_a.append(a )
        pred_e.append(e )
        pred_p.append(pl)
        pred_g.append(g )
    return pred_a, pred_e, pred_p, pred_g