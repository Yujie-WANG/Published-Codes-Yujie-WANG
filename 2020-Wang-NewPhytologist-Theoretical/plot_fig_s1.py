# import path and add new path to the list

from pylab import figure,rcParams

from photosynthesis import get_a_ci
from make_panel_label     import make_panel_labels

rcParams['font.family'] = 'serif'
rcParams['mathtext.fontset'] = 'dejavuserif'




# function to generate data for plotting
def get_e_a_dade(d_leaf=1.5, ca=40.0):
    e_leaf = 0.0
    e_de   = 1E-6
    list_e = []
    list_a = []
    list_s = []
    while True:
        g_sw  = e_leaf / (0.01*d_leaf)
        g_sc  = g_sw   / 1.6
        g_tar = g_sc   * 10.0
        ci,a  = get_a_ci(75.0, 120.0, 2.5, g_tar, ca, 25.0, 1000.0)
        g_sw_de  = (e_leaf+e_de) / (0.01*d_leaf)
        g_sc_de  = g_sw_de       / 1.6
        g_tar_de = g_sc_de       * 10.0
        ci,a_de  = get_a_ci(75.0, 120.0, 2.5, g_tar_de, ca, 25.0, 1000.0)
        dade     = (a_de - a) / e_de
        list_e.append(e_leaf*1E3)
        list_a.append(a         )
        list_s.append(dade*1E-3 )
        if e_leaf > 0.01:
            break
        e_leaf += 1E-4
    return list_e, list_a, list_s




le  ,la  ,ls   = get_e_a_dade()
le_d,la_d,ls_d = get_e_a_dade(d_leaf=3.0)
le_c,la_c,ls_c = get_e_a_dade(ca=80.0   )

fig = figure(figsize=(4,4), dpi=100)
fig.set_tight_layout(True)

ax1 = fig.add_subplot(1,1,1)


ax1.plot(le  ,ls  , "g:", label="$D$ = 1.5 kPa")
ax1.plot(le_d,ls_d, "g-", label="$D$ = 3.0 kPa")

ax1.set_xlim(0, 3)
ax1.set_ylim(0,16)

fig.text(0.02,0.19, "$\partial A/ \partial E$"     , color="green",
         fontsize=16, rotation=90, ha="right", va="bottom")
fig.text(0.02,0.96, "$\partial \Theta/ \partial E$", color="red"  ,
         fontsize=16, rotation=90, ha="right", va="top"   )
ax1.set_xlabel("$E_\mathrm{leaf}$ ($\mathrm{mmol\ m^{-2}\ s^{-1}}$)",
               color="black", fontsize=16)

fig.savefig("fig-s1-bak.svg", bbox_inches="tight")