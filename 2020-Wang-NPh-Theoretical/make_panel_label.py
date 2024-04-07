# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 08:25:18 2018

@author: Yujie
"""

from numpy import exp,log

csfont = {'fontname':'Serif', 'fontweight':'bold'}

# make labels
def make_panel_label(label,
                     ax,
                     loc="upper left",
                     paren=True,
                     logscale=False,
                     textbg="None",
                     fs=20,
                     distance_x=0.035,
                     distance_y=0.035):
    if paren:
        new_label = "(" + label + ")"
    if logscale==True:
        x_min = log( ax.get_xlim()[0] )
        x_max = log( ax.get_xlim()[1] )
        y_min = log( ax.get_ylim()[0] )
        y_max = log( ax.get_ylim()[1] )
    else:
        x_min = ax.get_xlim()[0]
        x_max = ax.get_xlim()[1]
        y_min = ax.get_ylim()[0]
        y_max = ax.get_ylim()[1]
    if( loc=="upper left" ):
        loc_x = x_min + (x_max-x_min) * distance_x
        loc_y = y_max - (y_max-y_min) * distance_y
        ha = "left"
        va = "top"
    elif( loc=="upper right" ):
        loc_x = x_max - (x_max-x_min) * distance_x
        loc_y = y_max - (y_max-y_min) * distance_y
        ha = "right"
        va = "top"
    elif( loc=="upper center" ):
        loc_x = x_min + (x_max-x_min) * 0.50
        loc_y = y_max - (y_max-y_min) * distance_y
        ha = "center"
        va = "top"
    elif( loc=="lower left" ):
        loc_x = x_min + (x_max-x_min) * distance_x
        loc_y = y_min + (y_max-y_min) * distance_y
        ha = "left"
        va = "bottom"
    if logscale==True:
        loc_x = exp(loc_x)
        loc_y = exp(loc_y)
    if textbg=="None":
        ax.text(loc_x, loc_y, new_label, fontsize=fs, ha=ha, va=va,
                **csfont)
    else:
        ax.text(loc_x, loc_y, new_label, fontsize=fs, ha=ha, va=va,
                **csfont,
                bbox=dict(boxstyle="square", facecolor=textbg, edgecolor="none"))

LIST_LABEL_UPPER = ["A", "B", "C", "D", "E", "F", "G",
                    "H", "I", "J", "K", "L", "M", "N",
                    "O", "P", "Q", "R", "S", "T", "U",
                    "V", "W", "X", "Y", "Z"]
LIST_LABEL_LOWER = ["a", "b", "c", "d", "e", "f", "g",
                    "h", "i", "j", "k", "l", "m", "n",
                    "o", "p", "q", "r", "s", "t", "u",
                    "v", "w", "x", "y", "z"]

def make_panel_labels(list_ax,
                      upper=True,
                      loc="upper left",
                      paren=True,
                      logscale=False,
                      textbg="None",
                      fs=20,
                      distance_x=0.035,
                      distance_y=0.035):
    if upper:
        for i in range(len(list_ax)):
            make_panel_label(LIST_LABEL_UPPER[i],
                             list_ax[i],
                             loc,
                             paren,
                             logscale,
                             textbg,
                             fs,
                             distance_x,
                             distance_y)
    else:
        for i in range(len(list_ax)):
            make_panel_label(LIST_LABEL_LOWER[i],
                             list_ax[i],
                             loc,
                             paren,
                             logscale,
                             textbg,
                             fs,
                             distance_x,
                             distance_y)