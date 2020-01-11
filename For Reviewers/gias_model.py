# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 09:30:30 2018

@author: yujie
"""

# import modules
from numpy import log,pi,sqrt

# class of an element
class diffusion_element:
    def __init__(self):
        self.l_hor_i = 0.1
        self.l_hor_o = 0.1
        self.l_ver_i = 1.0
        self.l_ver_o = 1.0
        self.y_hor_i = 1.0
        self.y_hor_o = 1.0
        self.y_ver_i = 1.0
        self.y_ver_o = 1.0
        self.p       = 0.0
    def load_values(self, values):
        self.l_hor_i = values[0]
        self.l_hor_o = values[1]
        self.l_ver_i = values[2]
        self.l_ver_o = values[3]
        self.y_hor_i = values[4]
        self.y_hor_o = values[5]
        self.y_ver_i = values[6]
        self.y_ver_o = values[7]

# class of a stoma model
class stoma_model():
    # define the class variables
    def __init__(self):
        self.r_ib = 100.0 # um
        self.r_st = 3.0   # um
        self.l_ia = 100.0 # um
        self.l_st = 20.0  # um
        self.l_bl = 400.0 # um
        self.mat_ia = []
        self.mat_bl = []
    # define the open in/out cylinder, applies to bl type
    def make_cylinder_bl(self,depth,printing=False):
        h_vert = depth * 0.05
        r_list = []
        a_list = []
        m_list = []
        h_list_i = [] # this is the ratio of r_2/r_1
        h_list_o = []
        h_diff_i = []
        h_diff_o = []
        v_list_i = []
        v_list_o = []
        # get the r_out, area, r_effective of the shells
        n_st = 1
        if self.r_st <= 5.0:
            r_list.append( self.r_st )
            a_list.append( pi * self.r_st ** 2.0 )
            m_list.append( sqrt(0.5) * self.r_st )
            n_st = 1
        else:
            r_list.append( self.r_st * 0.5 )
            a_list.append( pi * self.r_st ** 2.0 * 0.25 )
            m_list.append( sqrt(0.5) * self.r_st * 0.5 )
            r_list.append( self.r_st )
            a_list.append( pi * self.r_st ** 2.0 * 0.75 )
            m_list.append( sqrt(2.5) * self.r_st * 0.5)
            n_st = 2
        for i in range(1,11):
            r_list.append( (self.r_ib-self.r_st)*0.1*i + self.r_st )
            a_list.append( pi* (r_list[-1]**2.0 - r_list[-2]**2.0) )
            m_list.append( sqrt(0.5 * (r_list[-1]**2.0 + r_list[-2]**2.0)) )
        # get the horizontal and vertical distance of diffusion
        for i in range(len(r_list)):
            if i==0:
                h_list_i.append( 0.1 )
                h_list_o.append( m_list[i+1] / m_list[i] )
                h_diff_i.append( 0.0 )
                h_diff_o.append( 1.0 )
            elif i<len(r_list)-1:
                h_list_i.append( m_list[i] / m_list[i-1] )
                h_list_o.append( m_list[i+1] / m_list[i] )
                h_diff_i.append( 1.0 )
                h_diff_o.append( 1.0 )
            else:
                h_list_i.append( m_list[i] / m_list[i-1] )
                h_list_o.append( r_list[i] / m_list[i] )
                h_diff_i.append( 1.0 )
                h_diff_o.append( 0.0 )
        for i in range(20):
            if i==0:
                v_list_i.append( 0.5 * 0.05 * depth )
                v_list_o.append( 0.05 * depth )
            elif i<19:
                v_list_i.append( 0.05 * depth )
                v_list_o.append( 0.05 * depth )
            else:
                v_list_i.append( 0.05 * depth )
                v_list_o.append( 0.5 * 0.05 * depth )
        # make the matrix
        mat_ele = []
        ini_ele = []
        end_ele = []
        for i in range(len(r_list)+2):
            tmp_ele = diffusion_element()
            tmp_ele.p = 100.0
            ini_ele.append( tmp_ele )
        for i in range(len(r_list)+2):
            tmp_ele = diffusion_element()
            tmp_ele.p = 0.0
            end_ele.append( tmp_ele )
        mat_ele.append(ini_ele)
        for i in range(20):
            lay_ele = []
            tmp_ini = diffusion_element()
            tmp_end = diffusion_element()
            lay_ele.append( tmp_ini )
            if i==0:
                for j in range(len(r_list)):
                    tmp_ele = diffusion_element()
                    if j<n_st:
                        tmp_ele.load_values([h_list_i[j], h_list_o[j], v_list_i[i], v_list_o[i],
                                             h_diff_i[j], h_diff_o[j], 1.0        , 1.0        ])
                    else:
                        tmp_ele.load_values([h_list_i[j], h_list_o[j], v_list_i[i], v_list_o[i],
                                             h_diff_i[j], h_diff_o[j], 0.0        , 1.0        ])
                    lay_ele.append( tmp_ele )
            else:
                for j in range(len(r_list)):
                    tmp_ele = diffusion_element()
                    tmp_ele.load_values([h_list_i[j], h_list_o[j], v_list_i[i], v_list_o[i],
                                         h_diff_i[j], h_diff_o[j], 1.0        , 1.0        ])
                    lay_ele.append( tmp_ele )
            lay_ele.append( tmp_end )
            mat_ele.append( lay_ele )
        mat_ele.append(end_ele)
        # run the diffusion
        count = 0
        diff_in = 0.0
        diff_ou = 0.0
        j_in = 0.0
        g_in = 0.0
        while True:
            count += 1
            diff_tim = 0.5 / 100.0 * (depth * 0.05) ** 2.0
            diff_all = 0.0
            mat_dif = []
            # calculate the pressure increase
            for i in range(1,len(mat_ele)-1):
                lay_dif = []
                for j in range(1,len(mat_ele[0])-1):
                    diff  = 0.0
                    diff += (mat_ele[i-1][j].p - mat_ele[i][j].p) / mat_ele[i][j].l_ver_i * mat_ele[i][j].y_ver_i * a_list[j-1]
                    diff += (mat_ele[i+1][j].p - mat_ele[i][j].p) / mat_ele[i][j].l_ver_o * mat_ele[i][j].y_ver_o * a_list[j-1]
                    diff += 2.0 * pi * (mat_ele[i][j-1].p - mat_ele[i][j].p) / log(mat_ele[i][j].l_hor_i) * mat_ele[i][j].y_hor_i * h_vert
                    diff += 2.0 * pi * (mat_ele[i][j+1].p - mat_ele[i][j].p) / log(mat_ele[i][j].l_hor_o) * mat_ele[i][j].y_hor_o * h_vert
                    diff /= (a_list[j-1] * h_vert)
                    diff *= diff_tim
                    lay_dif.append( diff )
                    diff_all += abs(diff)
                mat_dif.append( lay_dif )
            # update the pressure
            for i in range(1,len(mat_ele)-1):
                for j in range(1,len(mat_ele[0])-1):
                    mat_ele[i][j].p += mat_dif[i-1][j-1]
            # calculate in and out
            diff_in = 0.0
            diff_ou = 0.0
            for i in range(1,len(mat_ele[0])-1):
                diff_in += (mat_ele[0][i].p - mat_ele[1][i].p) / (mat_ele[1][i].l_ver_i*1E-4) * mat_ele[1][i].y_ver_i * (a_list[i-1]*1E-8)
            for i in range(1,len(mat_ele[0])-1):
                diff_ou += (mat_ele[20][i].p - mat_ele[21][i].p) / (mat_ele[20][i].l_ver_o*1E-4) * mat_ele[20][i].y_ver_o * (a_list[i-1]*1E-8)
            diff_rat = abs(diff_in-diff_ou) / (abs(diff_in)+abs(diff_ou))
            if count==100 and printing==True:
                print(diff_in, diff_ou, diff_rat, "BL Type")
                count = 0
            if diff_rat < 1E-3 or count>100000:
                #print(diff_in, diff_ou, diff_rat, "BL Type")
                self.mat_bl = mat_ele
                print(count)
                break
        # diff in in unit of kPa / cm * cm2, j_in then cm2 /s * kPa *cm / (cm3 * kPa * mol / K * K) / cm2 = mol / cm2 / s
        j_in = 0.282 * diff_in / (8.314472*1E6*1E-3*298.13) / (pi * (self.r_ib*1E-4)**2.0)
        # g_in = j_in / d, d=100.0 kPa / 100.0 kPa, g_in in mol m-2 s-1 * density
        g_in = j_in * 1E4 * 1.013
        #print( "\nThe diff_in, diff_out and g_in are:", diff_in, diff_ou, g_in, "\n" )
        return g_in
    # define the open in/out cylinder, applies to ia type
    def make_cylinder_ia(self, depth, printing=False):
        h_vert = depth * 0.05
        r_list = []
        a_list = []
        m_list = []
        h_list_i = [] # this is the ratio of r_2/r_1
        h_list_o = []
        h_diff_i = []
        h_diff_o = []
        v_list_i = []
        v_list_o = []
        # get the r_out, area, r_effective of the shells
        n_st = 1
        if self.r_st <= 5.0:
            r_list.append( self.r_st )
            a_list.append( pi * self.r_st ** 2.0 )
            m_list.append( sqrt(0.5) * self.r_st )
            n_st = 1
        else:
            r_list.append( self.r_st * 0.5 )
            a_list.append( pi * self.r_st ** 2.0 * 0.25 )
            m_list.append( sqrt(0.5) * self.r_st * 0.5 )
            r_list.append( self.r_st )
            a_list.append( pi * self.r_st ** 2.0 * 0.75 )
            m_list.append( sqrt(2.5) * self.r_st * 0.5)
            n_st = 2
        for i in range(1,11):
            r_list.append( (self.r_ib-self.r_st)*0.1*i + self.r_st )
            a_list.append( pi* (r_list[-1]**2.0 - r_list[-2]**2.0) )
            m_list.append( sqrt(0.5 * (r_list[-1]**2.0 + r_list[-2]**2.0)) )
        # get the horizontal and vertical distance of diffusion
        for i in range(len(r_list)):
            if i==0:
                h_list_i.append( 0.1 )
                h_list_o.append( m_list[i+1] / m_list[i] )
                h_diff_i.append( 0.0 )
                h_diff_o.append( 1.0 )
            elif i<len(r_list)-1:
                h_list_i.append( m_list[i] / m_list[i-1] )
                h_list_o.append( m_list[i+1] / m_list[i] )
                h_diff_i.append( 1.0 )
                h_diff_o.append( 1.0 )
            else:
                h_list_i.append( m_list[i] / m_list[i-1] )
                h_list_o.append( r_list[i] / m_list[i] )
                h_diff_i.append( 1.0 )
                h_diff_o.append( 0.0 )
        for i in range(20):
            if i==0:
                v_list_i.append( 0.5 * 0.05 * depth )
                v_list_o.append( 0.05 * depth )
            elif i<19:
                v_list_i.append( 0.05 * depth )
                v_list_o.append( 0.05 * depth )
            else:
                v_list_i.append( 0.05 * depth )
                v_list_o.append( 0.5 * 0.05 * depth )
        # make the matrix
        mat_ele = []
        ini_ele = []
        end_ele = []
        for i in range(len(r_list)+2):
            tmp_ele = diffusion_element()
            if(i<=n_st):
                tmp_ele.p = 100.0
            else:
                tmp_ele.p = 0.0
            ini_ele.append( tmp_ele )
        for i in range(len(r_list)+2):
            tmp_ele = diffusion_element()
            tmp_ele.p = 0.0
            end_ele.append( tmp_ele )
        mat_ele.append(ini_ele)
        for i in range(20):
            lay_ele = []
            tmp_ini = diffusion_element()
            tmp_end = diffusion_element()
            lay_ele.append( tmp_ini )
            for j in range(len(r_list)):
                tmp_ele = diffusion_element()
                tmp_ele.load_values([h_list_i[j], h_list_o[j], v_list_i[i], v_list_o[i],
                                     h_diff_i[j], h_diff_o[j], 1.0        , 1.0        ])
                lay_ele.append( tmp_ele )
            lay_ele.append( tmp_end )
            mat_ele.append( lay_ele )
        mat_ele.append(end_ele)
        # run the diffusion
        count = 0
        while True:
            count += 1
            diff_tim = 0.5 / 100.0 * (depth * 0.05) ** 2.0
            diff_all = 0.0
            mat_dif = []
            # calculate the pressure increase
            for i in range(1,len(mat_ele)-1):
                lay_dif = []
                for j in range(1,len(mat_ele[0])-1):
                    diff  = 0.0
                    diff += (mat_ele[i-1][j].p - mat_ele[i][j].p) / mat_ele[i][j].l_ver_i * mat_ele[i][j].y_ver_i * a_list[j-1]
                    diff += (mat_ele[i+1][j].p - mat_ele[i][j].p) / mat_ele[i][j].l_ver_o * mat_ele[i][j].y_ver_o * a_list[j-1]
                    diff += 2.0 * pi * (mat_ele[i][j-1].p - mat_ele[i][j].p) / log(mat_ele[i][j].l_hor_i) * mat_ele[i][j].y_hor_i * h_vert
                    diff += 2.0 * pi * (mat_ele[i][j+1].p - mat_ele[i][j].p) / log(mat_ele[i][j].l_hor_o) * mat_ele[i][j].y_hor_o * h_vert
                    diff /= (a_list[j-1] * h_vert)
                    diff *= diff_tim
                    lay_dif.append( diff )
                    diff_all += abs(diff)
                mat_dif.append( lay_dif )
            # update the pressure
            for i in range(1,len(mat_ele)-1):
                for j in range(1,len(mat_ele[0])-1):
                    mat_ele[i][j].p += mat_dif[i-1][j-1]
            # calculate in and out
            diff_in = 0.0
            diff_ou = 0.0
            for i in range(1,n_st+1):
                diff_in += (mat_ele[0][i].p - mat_ele[1][i].p) / (mat_ele[1][i].l_ver_i*1E-4) * mat_ele[1][i].y_ver_i * (a_list[i-1]*1E-8)
            for i in range(n_st+1,len(mat_ele[0])-1):
                diff_ou += (mat_ele[1][i].p - mat_ele[0][i].p) / (mat_ele[1][i].l_ver_i*1E-4) * mat_ele[1][i].y_ver_i * (a_list[i-1]*1E-8)
            for i in range(1,len(mat_ele[0])-1):
                diff_ou += (mat_ele[20][i].p - mat_ele[21][i].p) / (mat_ele[20][i].l_ver_o*1E-4) * mat_ele[20][i].y_ver_o * (a_list[i-1]*1E-8)
            diff_rat = abs(diff_in-diff_ou) / (abs(diff_in)+abs(diff_ou))
            if count==100 and printing==True:
                print(diff_in, diff_ou, diff_rat, "IA Type")
                count = 0
            if diff_rat < 1E-3:
                #print(diff_in, diff_ou, diff_rat, "IA Type")
                self.mat_ia = mat_ele
                break
        # diff in in unit of kPa / cm * cm2, j_in then cm2 /s * kPa *cm / (cm3 * kPa * mol / K * K) / cm2 = mol / cm2 / s
        j_in = 0.282 * diff_in / (8.314472*1E6*1E-3*298.13) / (pi * (self.r_ib*1E-4)**2.0)
        # g_in = j_in / d, d=100.0 kPa / 100.0 kPa, g_in in mol m-2 s-1 * density
        g_in = j_in * 1E4 * 1.013
        #print( "\nThe diff_in, diff_out and g_in are:", diff_in, diff_ou, g_in, "\n" )
        return g_in