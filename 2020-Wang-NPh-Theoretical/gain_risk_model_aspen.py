# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 12:46:23 2019

@author: jesin
"""

from numpy import array,exp,log,zeros

from photosynthesis import get_a_ci,get_a_ci_seg




# gain-risk class
class gain_risk_model_aspen():
    def __init__(self):
        # soil
        self.p_ssat = 0.63*9.8*998.0*1E-6 # p_soil at saturation
        self.p_soil = self.p_ssat         # current p_soil
        self.c_ssat = 0.476               # swc at saturation
        self.c_curr = 0.476               # current swc
        self.b_ssat = 8.52                # unitless soil property
        self.k_ssat = 0.009               # m per hour
        self.k_rhiz = 5E8                 # rhizosphere conductance per basal
        # root
        self.k_root = 2228.0 # root conductance per basal area
        self.b_root = 1.15  # root weibull b
        self.c_root = 1.07  # root weibull c
        # stem
        self.k_stem = 4926.7 # stem conductance per basal area
        self.b_stem = 3.12  # stem weibull b
        self.c_stem = 2.64  # stem weibull c
        # leaf
        self.k_leaf = 5424.3 # leaf conductance per basal area
        self.b_leaf = 1.71  # stem weibull b
        self.c_leaf = 1.08  # stem weibull c
        self.k_sla  = 1.140  # leaf conductance per leaf area
        # tree level
        self.laba = 4758.5 # leaf area per basal area
        self.gaba = 1000.0 # crown area per basal area
        self.vmax = 61.74  # per leaf area
        self.jmax = 111.13 # per leaf area
        # height related
        self.h_soil = 0.0 # soil depth
        self.h_root = 0.0 # root depth
        self.h_stem = 1.0 # stem height
        self.h_leaf = 0.0 # leaf height
        # legacy
        self.l_root = zeros(20) # pressure history of 20 elements
        self.l_stem = zeros(20) # pressure history of 20 elements
        self.l_leaf = zeros(20) # pressure history of 20 elements
        # gias related
        self.c_cons = 0.0
        self.c_pows = 0.3
    
    def get_p_leaf(self, flow):
        # 1.1 assign the soil parameters
        p_ssat  = self.p_ssat
        c_ssat  = self.c_ssat
        b_ssat  = self.b_ssat
        k_rhiz  = self.k_rhiz
        p_soil  = self.p_soil
        tension = self.p_soil
        dp      = 0.0
        # 1.2 iterate through each shell to get dp
        for i in range(10):
            if tension>p_ssat:
                shell_t = c_ssat * (tension/p_ssat)**(-1.0/b_ssat)
            else:
                shell_t = c_ssat
            shell_f  = (shell_t/c_ssat) ** (2.0*b_ssat+3)
            shell_k  = k_rhiz * shell_f * log(10.0) / log((10.0-0.9*i)/(10.0-0.9*(i+1)))
            shell_k  = max(shell_k, 1E-12)
            dp      += flow / shell_k
            tension  = p_soil + dp
        p_rhiz  = tension
        tension = p_rhiz
        
        # 2.1 assign the root parameters
        b       = self.b_root
        c       = self.c_root
        k_root  = self.k_root
        h       = self.h_root
        legacy  = self.l_root
        dp      = 0.0
        # 2.2 iterate through each layer to get root dp including gravity
        for i in range(20):
            p = max(legacy[i],tension)
            f = exp( -1.0 * (p/b)**c )
            layer_k = k_root * 20.0 * f
            layer_k = max(layer_k, 1E-12)
            dp += flow / layer_k + 998.0*9.8*h*0.05*1E-6
            tension = p_rhiz + dp
        p_root  = tension
        
        # 3.1 assign the stem parameters
        b       = self.b_stem
        c       = self.c_stem
        k_stem  = self.k_stem
        h       = self.h_stem
        legacy  = self.l_stem
        dp      = 0.0
        # 3.2 iterate through each layer to get stem dp including gravity
        for i in range(20):
            p = max(legacy[i],tension)
            f = exp( -1.0 * (p/b)**c )
            layer_k = k_stem * 20.0 * f
            layer_k = max(layer_k, 1E-12)
            dp += flow / layer_k + 998.0*9.8*h*0.05*1E-6
            tension = p_root + dp
        p_stem  = tension
        
        # 4.1 assign the leaf parameters
        b       = self.b_leaf
        c       = self.c_leaf
        k_leaf  = self.k_leaf
        h       = self.h_leaf
        legacy  = self.l_leaf
        dp      = 0.0
        # 4.2 iterate through each layer to get stem dp including gravity
        for i in range(20):
            p = max(legacy[i],tension)
            f = exp( -1.0 * (p/b)**c )
            layer_k = k_leaf * 20.0 * f
            layer_k = max(layer_k, 1E-12)
            dp += flow / layer_k + 998.0*9.8*h*0.05*1E-6
            tension = p_stem + dp
        p_leaf  = tension
        
        return p_leaf
    
    def get_e_crit(self):
        p_crit = self.b_leaf * log(1000.0) ** (1.0/self.c_leaf)
        e_min  = 0.0
        e_max  = 100.0
        e_crit = 50.0
        while True:
            p = self.get_p_leaf(e_max)
            if p<p_crit:
                e_max *= 2.0
            else:
                break
        while True:
            e = 0.5 * (e_max+e_min)
            p = self.get_p_leaf(e)
            if abs(p-p_crit)<1E-3 or (e_max-e_min)<1E-3:
                e_crit = e
                break
            if p>p_crit:
                e_max  = e
            else:
                e_min  = e
        return e_crit
    
    def get_optima_dewar(self, d_leaf=1.5, ca=40.0, t_leaf=25.0, par=1000.0):
        e_crit = self.get_e_crit()
        p_crit = self.b_leaf * log(1000.0) ** (1.0/self.c_leaf)
        e_min  = 0.0
        e_max  = e_crit
        e_opt  = 0.0
        a_opt  = 0.0
        de     = 1.0
        while True:
            e   = 0.5 * (e_min+e_max)
            p   = self.get_p_leaf(e)
            f   = e * 0.0154321
            gh  = f / self.laba / d_leaf * 100.0
            gc  = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g   = gc * 10.0
            c,a = get_a_ci_seg(self.vmax,
                               self.jmax,
                               2.5,
                               g,
                               ca,
                               t_leaf,
                               par)
            e_de = e + de
            p_de = self.get_p_leaf(e_de)
            f_de = e_de * 0.0154321
            gh   = f_de / self.laba / d_leaf * 100.0
            gc   = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g_de = gc * 10.0
            c_de,a_de = get_a_ci_seg(self.vmax,
                                     self.jmax,
                                     2.5,
                                     g_de,
                                     ca,
                                     t_leaf,
                                     par)
            optimizer = a_de * (p_crit-p_de) - a * (p_crit-p)
            if (e_max-e_min)<1.0:
                e_opt = e
                p_opt = self.get_p_leaf(e)
                a_opt = a * p_opt / p_crit
                break
            if optimizer>0:
                e_min = e
            else:
                e_max = e
        return a_opt,e_opt,p_opt
    
    def get_optima_dewar_mod(self, d_leaf=1.5, ca=40.0, t_leaf=25.0, par=1000.0):
        e_crit = self.get_e_crit()
        p_crit = self.b_leaf * log(1000.0) ** (1.0/self.c_leaf)
        e_min  = 0.0
        e_max  = e_crit
        e_opt  = 0.0
        a_opt  = 0.0
        de     = 1.0
        while True:
            e   = 0.5 * (e_min+e_max)
            p   = self.get_p_leaf(e)
            f   = e * 0.0154321
            gh  = f / self.laba / d_leaf * 100.0
            gc  = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g   = gc * 10.0
            c,a = get_a_ci_seg(self.vmax,
                               self.jmax,
                               2.5,
                               g,
                               ca,
                               t_leaf,
                               par)
            e_de = e + de
            p_de = self.get_p_leaf(e_de)
            f_de = e_de * 0.0154321
            gh   = f_de / self.laba / d_leaf * 100.0
            gc   = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g_de = gc * 10.0
            c_de,a_de = get_a_ci_seg(self.vmax,
                                     self.jmax,
                                     2.5,
                                     g_de,
                                     ca,
                                     t_leaf,
                                     par)
            optimizer = a_de * (p_crit-p_de) - a * (p_crit-p)
            if (e_max-e_min)<1.0:
                e_opt = e
                p_opt = self.get_p_leaf(e)
                a_opt = a
                break
            if optimizer>0:
                e_min = e
            else:
                e_max = e
        return a_opt,e_opt,p_opt
    
    def get_optima_eller(self, d_leaf=1.5, ca=40.0, t_leaf=25.0, par=1000.0):
        e_crit = self.get_e_crit()
        e_min  = 0.0
        e_max  = e_crit
        e_opt  = 0.0
        a_opt  = 0.0
        de     = 1.0
        while True:
            e   = 0.5 * (e_min+e_max)
            p   = self.get_p_leaf(e)
            f   = e * 0.0154321
            gh  = f / self.laba / d_leaf * 100.0
            gc  = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g   = gc * 10.0
            c,a = get_a_ci_seg(self.vmax,
                               self.jmax,
                               2.5,
                               g,
                               ca,
                               t_leaf,
                               par)
            e_de = e + de
            p_de = self.get_p_leaf(e_de)
            f_de = e_de * 0.0154321
            gh   = f_de / self.laba / d_leaf * 100.0
            gc   = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g_de = gc * 10.0
            c_de,a_de = get_a_ci_seg(self.vmax,
                                     self.jmax,
                                     2.5,
                                     g_de,
                                     ca,
                                     t_leaf,
                                     par)
            e_df = e_de + de
            p_df = self.get_p_leaf(e_df)
            m    = de / (p_de - p   )
            m_de = de / (p_df - p_de)
            optimizer = a_de * m_de - a * m
            if (e_max-e_min)<1.0:
                e_opt = e
                a_opt = a
                p_opt = self.get_p_leaf(e)
                break
            if optimizer>0:
                e_min = e
            else:
                e_max = e
        return a_opt,e_opt,p_opt
    
    def get_optima_lambda(self, lambd=0.0, d_leaf=1.5, ca=40.0, t_leaf=25.0, par=1000.0):
        e_crit = self.get_e_crit()
        e_min  = 0.0
        e_max  = e_crit
        e_opt  = 0.0
        a_opt  = 0.0
        de     = 1.0
        while True:
            e   = 0.5 * (e_min+e_max)
            f   = e * 0.0154321
            gh  = f / self.laba / d_leaf * 100.0
            gc  = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g   = gc * 10.0
            c,a = get_a_ci_seg(self.vmax,
                               self.jmax,
                               2.5,
                               g,
                               ca,
                               t_leaf,
                               par)
            e_de = e + de
            f_de = e_de * 0.0154321
            gh   = f_de / self.laba / d_leaf * 100.0
            gc   = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g_de = gc * 10.0
            c_de,a_de = get_a_ci_seg(self.vmax,
                                     self.jmax,
                                     2.5,
                                     g_de,
                                     ca,
                                     t_leaf,
                                     par)
            optimizer = (a_de - a) / de
            if (e_max-e_min)<1.0:
                e_opt = e
                a_opt = a
                p_opt = self.get_p_leaf(e)
                break
            if optimizer>lambd:
                e_min = e
            else:
                e_max = e
        return a_opt,e_opt,p_opt
    
    def get_optima_prentice(self, ce=1.0, cv=1.0, d_leaf=1.5, ca=40.0, t_leaf=25.0, par=1000.0):
        e_crit = self.get_e_crit()
        e_min  = 0.0
        e_max  = e_crit
        e_opt  = 0.0
        a_opt  = 0.0
        de     = 1.0
        while True:
            e   = 0.5 * (e_min+e_max)
            f   = e * 0.0154321
            gh  = f / self.laba / d_leaf * 100.0
            gc  = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g   = gc * 10.0
            c,a = get_a_ci_seg(self.vmax,
                               self.jmax,
                               2.5,
                               g,
                               ca,
                               t_leaf,
                               par)
            e_de = e + de
            f_de = e_de * 0.0154321
            gh   = f_de / self.laba / d_leaf * 100.0
            gc   = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g_de = gc * 10.0
            c_de,a_de = get_a_ci_seg(self.vmax,
                                     self.jmax,
                                     2.5,
                                     g_de,
                                     ca,
                                     t_leaf,
                                     par)
            optimizer = a_de / (ce*e_de + cv*self.vmax) - a / (ce*e + cv*self.vmax)
            if (e_max-e_min)<1.0:
                e_opt = e
                a_opt = a
                p_opt = self.get_p_leaf(e)
                break
            if optimizer>0:
                e_min = e
            else:
                e_max = e
        return a_opt,e_opt,p_opt
    
    def get_optima_sperry(self, d_leaf=1.5, ca=40.0, t_leaf=25.0, par=1000.0):
        # 1. calculate the p_crit @ layer_f = 1E-6
        e_crit = self.get_e_crit()
        de = 1.0
        # 2. increase the e stepwise
        list_e = []
        list_k = []
        list_a = []
        list_p = []
        for i in range(101):
            e   = i * 0.01 * e_crit
            p   = self.get_p_leaf(e)
            f   = e * 0.0154321
            gh  = f / self.laba / d_leaf * 100.0
            gc  = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g   = gc * 10.0
            c,a = get_a_ci(self.vmax,
                           self.jmax,
                           2.5,
                           g,
                           ca,
                           t_leaf,
                           par)
            e_de = e + de
            p_de = self.get_p_leaf(e_de)
            k = de / (p_de-p)
            list_e.append(e)
            list_k.append(k)
            list_a.append(a)
            list_p.append(p)
        # 3. extend the lists
        gain = array(list_a)/max(list_a)
        risk = 1.0 - array(list_k)/max(list_k)
        prof = list( gain - risk )
        opt_site = prof.index(max(prof))
        opt_a = list_a[opt_site]
        opt_e = list_e[opt_site]
        opt_p = list_p[opt_site]
        return opt_a, opt_e, opt_p
    
    def get_optima_wang(self, d_leaf=1.5, ca=40.0, t_leaf=25.0, par=1000.0):
        e_crit = self.get_e_crit()
        e_min  = 0.0
        e_max  = e_crit
        e_opt  = 0.0
        a_opt  = 0.0
        de     = 1.0
        while True:
            e   = 0.5 * (e_min+e_max)
            f   = e * 0.0154321
            gh  = f / self.laba / d_leaf * 100.0
            gc  = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g   = gc * 10.0
            c,a = get_a_ci_seg(self.vmax,
                               self.jmax,
                               2.5,
                               g,
                               ca,
                               t_leaf,
                               par)
            e_de = e + de
            f_de = e_de * 0.0154321
            gh   = f_de / self.laba / d_leaf * 100.0
            gc   = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g_de = gc * 10.0
            c_de,a_de = get_a_ci_seg(self.vmax,
                                     self.jmax,
                                     2.5,
                                     g_de,
                                     ca,
                                     t_leaf,
                                     par)
            optimizer = a_de * (e_crit-e_de) - a * (e_crit-e)
            if (e_max-e_min)<1.0:
                e_opt = e
                a_opt = a
                p_opt = self.get_p_leaf(e)
                break
            if optimizer>0:
                e_min = e
            else:
                e_max = e
        return a_opt,e_opt,p_opt
    
    def get_optima_wap(self, aa=0.1, bb=0.1, d_leaf=1.5, ca=40.0, t_leaf=25.0, par=1000.0):
        e_crit = self.get_e_crit()
        e_min  = 0.0
        e_max  = e_crit
        e_opt  = 0.0
        a_opt  = 0.0
        de     = 1.0
        while True:
            e   = 0.5 * (e_min+e_max)
            p   = self.get_p_leaf(e)
            f   = e * 0.0154321
            gh  = f / self.laba / d_leaf * 100.0
            gc  = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g   = gc * 10.0
            c,a = get_a_ci_seg(self.vmax,
                               self.jmax,
                               2.5,
                               g,
                               ca,
                               t_leaf,
                               par)
            e_de = e + de
            p_de = self.get_p_leaf(e_de)
            f_de = e_de * 0.0154321
            gh   = f_de / self.laba / d_leaf * 100.0
            gc   = gh / 1.6 / (1.0 + self.c_cons * gh**self.c_pows)
            g_de = gc * 10.0
            c_de,a_de = get_a_ci_seg(self.vmax,
                                     self.jmax,
                                     2.5,
                                     g_de,
                                     ca,
                                     t_leaf,
                                     par)
            optimizer = a_de - aa*p_de**2.0 - bb*p_de - a + aa*p**2.0 + bb*p
            if (e_max-e_min)<1.0:
                e_opt = e
                a_opt = a
                p_opt = self.get_p_leaf(e)
                break
            if optimizer>0:
                e_min = e
            else:
                e_max = e
        return a_opt,e_opt,p_opt
