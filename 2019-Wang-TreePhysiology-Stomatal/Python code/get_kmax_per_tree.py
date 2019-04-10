# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 11:52:13 2018

@author: Yujie Wang
"""

from numpy  import exp
from pandas import read_csv

def get_optimal_kmax(list_ppd, list_pmd, list_emd):
    # assign the VCs
    root_b = 1.879
    root_c = 2.396
    stem_b = 2.238
    stem_c = 9.380
    leaf_b = 1.897
    leaf_c = 2.203
    # define the kmax to optimize
    maxk = 1E6
    mink = 1E-9
    while True:
        if maxk-mink<0.2:
            break
        kmax = 0.5 * (maxk+mink)
        root_kmax = kmax * 1.863013458
        stem_kmax = kmax * 4.11942096
        leaf_kmax = kmax * 4.535504046
        diff_p    = 0.0
        for i in range(len(list_ppd)):
            ppd = list_ppd[i]
            pmd = list_pmd[i]
            emd = list_emd[i]
            tmp = ppd
            for j in range(50):
                k = max(1E-12, exp(-(tmp/root_b)**root_c) * root_kmax * 50.0)
                tmp += emd / k
            for j in range(50):
                k = max(1E-12, exp(-(tmp/stem_b)**stem_c) * stem_kmax * 50.0)
                tmp += emd / k
            for j in range(50):
                k = max(1E-12, exp(-(tmp/leaf_b)**leaf_c) * leaf_kmax * 50.0)
                tmp += emd / k
            diff_p += (tmp-pmd)**2.0
        kmax += 0.1
        root_kmax = kmax * 1.863013458
        stem_kmax = kmax * 4.11942096
        leaf_kmax = kmax * 4.535504046
        diff_pdk  = 0.0
        for i in range(len(list_ppd)):
            ppd = list_ppd[i]
            pmd = list_pmd[i]
            emd = list_emd[i]
            tmp = ppd
            for j in range(50):
                k = max(1E-12, exp(-(tmp/root_b)**root_c) * root_kmax * 50.0)
                tmp += emd / k
            for j in range(50):
                k = max(1E-12, exp(-(tmp/stem_b)**stem_c) * stem_kmax * 50.0)
                tmp += emd / k
            for j in range(50):
                k = max(1E-12, exp(-(tmp/leaf_b)**leaf_c) * leaf_kmax * 50.0)
                tmp += emd / k
            diff_pdk += (tmp-pmd)**2.0
        #print(diff_p, "\t", diff_pdk, "\t", kmax)
        if diff_pdk-diff_p>0:
            maxk = kmax-0.1
        else:
            mink = kmax-0.1
    return kmax

data = read_csv("ppe.txt", delimiter="\t")

for no in range(1,11):
    # read the data subset
    data_sub = data.query("TreeNo==%d" % no)
    # interate through the subset
    list_ppd = data_sub["Ppd"  ].values
    list_pmd = data_sub["Pmd"  ].values
    list_emd = data_sub["Etree"].values
    kmax = get_optimal_kmax(list_ppd, list_pmd, list_emd)
    print("%d\t%.3f" % (no,kmax))