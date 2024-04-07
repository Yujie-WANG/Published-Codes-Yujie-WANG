# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 12:46:23 2019

@author: jesin
"""

from numpy import array,exp,log

from photosynthesis import get_a_ci,get_a_ci_seg




# gain-risk class
class gain_risk_model_other():
    def __init__(self):
        # tree no soil type given
        self.k_tree = 10.0 # tree conductance per basal area
        self.b_tree = 2.0   # tree weibull b
        self.c_tree = 5.0   # tree weibull c
        # tree level
        self.laba = 1000.0 # leaf area per basal area
        self.vmax = 61.74  # per leaf area
        self.jmax = 111.13 # per leaf area
        # gias related
        self.c_cons = 0.0
        self.c_pows = 0.3
    
    def get_p_leaf(self, flow):
        p_soil  = self.p_soil
        tension = p_soil
        
        # 2.1 assign the tree parameters
        b       = self.b_tree
        c       = self.c_tree
        k_tree  = self.k_tree
        dp      = 0.0
        # 2.2 iterate through each layer to get tree dp including gravity
        for i in range(20):
            p = tension
            f = exp( -1.0 * (p/b)**c )
            layer_k = k_tree * 20.0 * f
            layer_k = max(layer_k, 1E-12)
            dp += flow / layer_k
            tension = p_soil + dp
        p_leaf  = tension
        
        return p_leaf
    
    def get_e_crit(self):
        p_crit = self.b_tree * log(1000.0) ** (1.0/self.c_tree)
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
        p_crit = self.b_tree * log(1000.0) ** (1.0/self.c_tree)
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
                g_opt = g * 1.6
                break
            if optimizer>0:
                e_min = e
            else:
                e_max = e
        return a_opt,g_opt,p_opt
    
    def get_optima_dewar_mod(self, d_leaf=1.5, ca=40.0, t_leaf=25.0, par=1000.0):
        e_crit = self.get_e_crit()
        p_crit = self.b_tree * log(1000.0) ** (1.0/self.c_tree)
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
                g_opt = g * 1.6
                break
            if optimizer>0:
                e_min = e
            else:
                e_max = e
        return a_opt,g_opt,p_opt
    
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
                g_opt = g * 1.6
                break
            if optimizer>0:
                e_min = e
            else:
                e_max = e
        return a_opt,g_opt,p_opt
    
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
                g_opt = g * 1.6
                p_opt = self.get_p_leaf(e)
                break
            if optimizer>lambd:
                e_min = e
            else:
                e_max = e
        return a_opt,g_opt,p_opt
    
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
                g_opt = g * 1.6
                p_opt = self.get_p_leaf(e)
                break
            if optimizer>0:
                e_min = e
            else:
                e_max = e
        return a_opt,g_opt,p_opt
    
    def get_optima_sperry(self, d_leaf=1.5, ca=40.0, t_leaf=25.0, par=1000.0):
        # 1. calculate the p_crit @ layer_f = 1E-6
        e_crit = self.get_e_crit()
        de = 1.0
        # 2. increase the e stepwise
        list_e = []
        list_k = []
        list_a = []
        list_p = []
        list_g = []
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
            list_g.append(g * 1.6)
        # 3. extend the lists
        gain = array(list_a)/max(list_a)
        risk = 1.0 - array(list_k)/max(list_k)
        prof = list( gain - risk )
        opt_site = prof.index(max(prof))
        opt_a = list_a[opt_site]
        opt_e = list_e[opt_site]
        opt_p = list_p[opt_site]
        opt_g = list_g[opt_site]
        return opt_a, opt_g, opt_p
    
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
                g_opt = g * 1.6
                p_opt = self.get_p_leaf(e)
                break
            if optimizer>0:
                e_min = e
            else:
                e_max = e
        return a_opt,g_opt,p_opt
