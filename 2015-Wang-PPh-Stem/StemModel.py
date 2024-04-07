from numpy import *
from pylab import *

Alpha = 0.05
F_BUB = Alpha
F_WAT = 1.0 - Alpha
P_INI = 0.3
P_OUT = 1.0

VOLUM_FACTOR = 1.0
TORCU_FACTOR = 1.0
ALPHA_FACTOR = 1.0
BLOCK_FACTOR = VOLUM_FACTOR * TORCU_FACTOR * ALPHA_FACTOR

R_O2 = 0.20
R_N2 = 0.80

#HENRY"S LAW
KCP_O = 1.30E-6 #mol.cm-3.atm-1
KCP_N = 6.10E-7
KCP_A = R_N2*KCP_N + R_O2*KCP_O
RT = 298.0 * 82.1 # atm.cm3.mol-1

VESSEL_FACTOR = ((1.0-Alpha)*KCP_A) / (Alpha/RT + (1.0-Alpha)*KCP_A)
print VESSEL_FACTOR

# FICK"S LAW
#D_O2 = 2.10E-5 # cm2.s-1
#D_N2 = 1.88E-5 # cm2.s-1
D_O2 = 1.00E-6 # cm2.s-1
D_N2 = 1.00E-6 # cm2.s-1
D_ALL = D_O2*R_O2 + D_N2*R_N2

#FICK"S LAW
R_IN  = 0.00 #cm
R_OUT = 0.40 #cm
R_BARK = 0.45 #cm
#print R_OUT - sqrt((R_IN**2 + R_OUT**2)/2.0)
#print pi * (R_OUT**2 - R_IN**2)

N_SHELL = 20
T_INTERVAL = 1 #s

L_SHELL = (R_OUT-R_IN) / N_SHELL

R_SHELL_OUT = []
R_SHELL_IN  = []
P_SHELL = []
for i in range(N_SHELL):
    R_SHELL_OUT.append(R_OUT- i*L_SHELL)
    R_SHELL_IN.append(R_OUT - (i+1)*L_SHELL)
    P_SHELL.append(P_INI)
S_SHELL = []
M_SHELL = []
for i in range(N_SHELL):
    S_SHELL.append(pi * (R_SHELL_OUT[i]**2 - R_SHELL_IN[i]**2) * VOLUM_FACTOR)
    M_SHELL.append(sqrt((R_SHELL_OUT[i]**2 + R_SHELL_IN[i]**2)/2))
M_SHELL_MM = array(M_SHELL) * 10.0
SUM_S_SHELL = sum(S_SHELL)
SB3_S_SHELL = sum(S_SHELL[0:3])
SE3_S_SHELL = sum(S_SHELL[-3:])
    
K_SHELL = []
for i in range(N_SHELL):
    if(i == 0):
        tmp_K = 2.0*pi*D_ALL/log(R_OUT/M_SHELL[0]) * BLOCK_FACTOR * T_INTERVAL
    else:
        tmp_K = 2.0*pi*D_ALL/log(M_SHELL[i-1]/M_SHELL[i]) * BLOCK_FACTOR * T_INTERVAL
    K_SHELL.append(tmp_K)

# MAIN FUNCTION
RESULTS = []
TIME = 0
"""
R_FILE = open("./Stem_Diffusion.txt","w+")
tmpstr = "Time & Shells"
for i in M_SHELL:
    tmpstr += "\t"
    tmpstr += str(i)
tmpstr += "\r\n"
R_FILE.write(tmpstr)
tmpstr = str(TIME)
for i in P_SHELL:
    tmpstr += "\t"
    tmpstr += str(i)
tmpstr += "\r\n"
R_FILE.write(tmpstr)
"""
ALL_IN = 0.0
ALL_LAST = 0.0
TIME_LIMIT = 5*24*60*60
while(TIME <= TIME_LIMIT):
    MOLES_SHELL = []
    TIME += T_INTERVAL
    print TIME_LIMIT - TIME + 1
    for i in range(N_SHELL):
        if(i == 0):
            IO = K_SHELL[i]*(P_OUT-P_SHELL[i])*KCP_A
            ALL_IN += IO
            MOLES_SHELL.append(IO)
        else:
            IO = K_SHELL[i]*(P_SHELL[i-1]-P_SHELL[i])*KCP_A
            MOLES_SHELL.append(IO)
    for i in range(N_SHELL):
        if(i == N_SHELL-1):
            tmpIO = MOLES_SHELL[i]
            P_SHELL[i] += tmpIO/(S_SHELL[i]*(1.0-Alpha)*VOLUM_FACTOR) / KCP_A * VESSEL_FACTOR
        else:
            tmpIO = MOLES_SHELL[i] - MOLES_SHELL[i+1]
            P_SHELL[i] += tmpIO/(S_SHELL[i]*(1.0-Alpha)*VOLUM_FACTOR) / KCP_A * VESSEL_FACTOR
    sum_pressure = 0.0
    sb3_pressure = 0.0
    se3_pressure = 0.0
    for i in range(N_SHELL):
        if(i <= 2):
            sb3_pressure += P_SHELL[i]*S_SHELL[i]
        if(i >= N_SHELL-3):
            se3_pressure += P_SHELL[i]*S_SHELL[i]
        sum_pressure += P_SHELL[i]*S_SHELL[i]
    average_BP = sum_pressure / SUM_S_SHELL
    average_B3 = sb3_pressure / SB3_S_SHELL
    average_E3 = se3_pressure / SE3_S_SHELL
    RESULTS.append([TIME,P_SHELL[:],[average_BP,average_B3,average_E3]])
all_in = 0.0
for i in range(N_SHELL):
    tmp = S_SHELL[i] * (P_SHELL[i] - P_INI) * (Alpha/RT + (1.0-Alpha)*KCP_A)
    all_in += tmp
print ALL_IN,all_in
"""
    tmpstr = str(TIME)
    for i in P_SHELL:
        tmpstr += "\t"
        tmpstr += str(i)
    tmpstr += "\r\n"
    R_FILE.write(tmpstr)
#R_FILE.close()
"""

# Time versus average BUBBLE PRESSURE

# Choose data to plot
Y1 = array(RESULTS[60*60 - 1][1]) * 101.3
Y2 = array(RESULTS[3*60*60 - 1][1]) * 101.3
Y3 = array(RESULTS[8*60*60 - 1][1]) * 101.3
Y4 = array(RESULTS[24*60*60 - 1][1]) * 101.3

figure(figsize=(12,12),dpi=100)

subplot(2,2,1)
plot(M_SHELL_MM,Y1,color="black",linewidth=1.5)
ylim(0,105)
title("Pressure Distribution after 1 h")
#xlabel("Distance to center of the stem (cm)")
ylabel("Bubble pressure in shells(kPa)",fontsize=16)
text(0.005,103,"A",fontsize=16,fontweight="bold",horizontalalignment="left",verticalalignment="top")
ax = gca()
for axis in ax.spines:
    ax.spines[axis].set_linewidth(1.5)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontsize(14)

subplot(2,2,2)
plot(M_SHELL_MM,Y2,color="black",linewidth=1.5)
ylim(0,105)
title("Pressure Distribution after 3 h")
#xlabel("Distance to center of the stem (cm)")
#ylabel("Bubble pressure in embolized vessels(kPa)")
text(0.005,103,"B",fontsize=16,fontweight="bold",horizontalalignment="left",verticalalignment="top")
ax = gca()
for axis in ax.spines:
    ax.spines[axis].set_linewidth(1.5)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontsize(14)

subplot(2,2,3)
plot(M_SHELL_MM,Y3,color="black",linewidth=1.5)
ylim(0,105)
title("Pressure Distribution after 8 h")
xlabel("Distance to center of the stem (mm)",fontsize=16)
ylabel("Bubble pressure in shells(kPa)",fontsize=16)
text(0.005,103,"C",fontsize=16,fontweight="bold",horizontalalignment="left",verticalalignment="top")
ax = gca()
for axis in ax.spines:
    ax.spines[axis].set_linewidth(1.5)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontsize(14)

subplot(2,2,4)
plot(M_SHELL_MM,Y4,color="black",linewidth=1.5)
ylim(0,105)
title("Pressure Distribution after 24 h")
xlabel("Distance to center of the stem (mm)",fontsize=16)
#ylabel("Bubble pressure in embolized vessels(kPa)")
text(0.005,103,"D",fontsize=16,fontweight="bold",horizontalalignment="left",verticalalignment="top")
ax = gca()
for axis in ax.spines:
    ax.spines[axis].set_linewidth(1.5)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontsize(14)

savefig("./Stem-Model-Time-Examples.pdf",dpi=100)

figure(figsize=(6,10),dpi=100)

subplot(2,1,1)
XX = []
YY = []
YB = []
YE = []
for i in RESULTS:
    XX.append(i[0]/60.0/60.0)
    YY.append(i[2][0]*101.3)
    YB.append(i[2][1]*101.3)
    YE.append(i[2][2]*101.3)
plot(XX,YY,color="black",linewidth=1.5,linestyle="-" ,label="All")
plot(XX,YB,color="black",linewidth=1.5,linestyle="--",label="First 3 shells")
plot(XX,YE,color="black",linewidth=1.5,linestyle="-.",label="Last 3 shells")
xlabel("Time (h)",fontsize=16)
ylabel("Average bubble pressure (kPa)",fontsize=16)
xlim(0,120)
text(0.5,109,"A",fontsize=16,fontweight="bold",horizontalalignment="left",verticalalignment="top")
legend(loc="lower right",frameon=False)
ax = gca()
for axis in ax.spines:
    ax.spines[axis].set_linewidth(1.5)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontsize(14)

subplot(2,1,2)
SX = []
SP = []
SY = []
for i in range(len(XX)-1):
    tmp_k = (YY[i+1]-YY[i]) / (XX[i+1]-XX[i])
    tmp_t = (101.3-YY[i]) / tmp_k
    SX.append(XX[i])
    SP.append(YY[i])
    SY.append(tmp_t)
print len(SY),len(SP),len(SX)
plot(SX,SP,color="black",linewidth=1.5,linestyle="-" ,label="Bubble Pressure")
xlabel("Time (h)",fontsize=16)
ylabel("Average Bubble Pressure (kPa)",fontsize=16)
text(0.5,109,"B",fontsize=16,fontweight="bold",horizontalalignment="left",verticalalignment="top")
legend(loc="lower left",frameon=False)
xlim(0,60)
ax = gca()
for axis in ax.spines:
    ax.spines[axis].set_linewidth(1.5)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontsize(14)
twinx()
plot(SX,SY,color="black",linewidth=1.5,linestyle="--",label="Time Constant")
ylabel(r"Time Constant (h$^{-1}$)",fontsize=16)
legend(loc="lower right",frameon=False)
xlim(0,60)
ax = gca()
for axis in ax.spines:
    ax.spines[axis].set_linewidth(1.5)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontsize(14)
savefig("./Stem-Time-PB.pdf",dpi=100)
