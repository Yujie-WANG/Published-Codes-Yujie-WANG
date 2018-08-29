# This file runs under python3
from numpy import array,linspace
from pylab import figure

# make labels
def make_panel_label(label, ax, loc="upper left", fs=20):
    if( loc=="upper left" ):
        loc_x = ax.get_xlim()[0] + (ax.get_xlim()[1]-ax.get_xlim()[0]) * 0.02
        loc_y = ax.get_ylim()[1] - (ax.get_ylim()[1]-ax.get_ylim()[0]) * 0.02
    ax.text(loc_x, loc_y, label, fontsize=fs, ha="left", va="top")

# make analysis
BUBBLE_PRESSURE = linspace(20, 100, 5,endpoint=True)
WATER_PRESSURET = linspace(0,110,1101,endpoint=True)
CAPILARY_PRESSURE = 7.1

KRATIO_LIST_PB = []
for PLC in [0.6,]:#PLC_LIST:
    for PB in BUBBLE_PRESSURE:
        klist = []
        for PW in WATER_PRESSURET:
            if((CAPILARY_PRESSURE+PW) > PB):
                LW = 1.0 - PB/(CAPILARY_PRESSURE+PW)
            else:
                LW = 0.0
            K = 2*LW / (1+LW**2)
            klist.append(K*PLC + 1.0 - PLC)
        K_ARRAY = array(klist)
        KRATIO_LIST_PB.append(K_ARRAY)

# analysis of fact
WATER_PRESSUREF = linspace(10.0,101,901,endpoint=True)
PLC = 0.6
PB = 50.0
KE = 0.0002
klist = []
hlist = []
elist = []
for PW in WATER_PRESSUREF:
	if((CAPILARY_PRESSURE+PW) > PB):
		LW = 1.0 - PB/(CAPILARY_PRESSURE+PW)
	else:
		LW = 0.0
	K = 2*LW / (1+LW**2)
	E = KE * (101.3 - CAPILARY_PRESSURE - PW) * 101.3/(CAPILARY_PRESSURE+PW)
	klist.append(K*PLC + 1.0 - PLC)
	hlist.append(K*PLC*0.2 + 1.0 - PLC)
	elist.append(E)
	#print(K,"\t",E)
K_ARRAY = array(klist)
H_ARRAY = array(hlist)
E_ARRAY = array(elist)
F_ARRAY = H_ARRAY - E_ARRAY

colors = ["grey", "black","black","black", "black"]
linestyles = ["-", "--", ":", "-.", "-"]

fig = figure(figsize=(10,5), dpi=100)
fig.set_tight_layout(True)

ax1 = fig.add_subplot(1,2,1)
for i in range(len(KRATIO_LIST_PB)):
    ax1.plot(WATER_PRESSURET, KRATIO_LIST_PB[i], linestyle=linestyles[i], color=colors[i], label=str("$P_\mathrm{B}$ = "+str(BUBBLE_PRESSURE[i])))
#ax1.set_xlabel("Absolute Pressure (kPa)", fontsize=18)
ax1.set_xlabel("Absolute Water Pressure (kPa)", fontsize=18)
ax1.set_ylabel("$k_\mathrm{h}$ Ratio", fontsize=18)
ax1.set_xlim(0, 120)
ax1.set_ylim(0, 1  )
ax1.legend(loc="lower right", ncol=2)

ax2 = fig.add_subplot(1,2,2)
ax2.plot(WATER_PRESSUREF, K_ARRAY, color="black", linestyle="--", label="Theory"     )
ax2.plot(WATER_PRESSUREF, E_ARRAY, color="black", linestyle=":" , label="Air Flow"   )
ax2.plot(WATER_PRESSUREF, F_ARRAY, color="black", linestyle="-" , label="Observation")
ax2.set_xlabel("Absolute Water Pressure (kPa)", fontsize=18)
#ax2.set_ylabel("$k_\mathrm{h}$ Ratio", fontsize=18)
ax2.set_xlim(0, 120)
ax2.set_ylim(0, 1.0)
ax2.legend(bbox_to_anchor=(0.62, 0.05), loc="lower left")

make_panel_label("(a)", ax1);
make_panel_label("(b)", ax2);

fig.savefig("theory-fact.eps", bbox_inches="tight")
fig.savefig("theory-fact.jpg", bbox_inches="tight")
fig.savefig("theory-fact.pdf", bbox_inches="tight")