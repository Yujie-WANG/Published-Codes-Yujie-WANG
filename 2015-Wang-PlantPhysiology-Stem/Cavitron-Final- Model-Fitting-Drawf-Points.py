import numpy
from scipy.optimize import  leastsq

# Central PLC, basal end parameter, distal end parameter, bubble pressure and suraface tension
NEWR = []
OLDR = []
center_PLC = 50
a_para = 0
b_para = 0
BP = 50

KCP_O = 1.30E-3 #mol.L-3.atm-1
KCP_N = 6.10E-4
KCP_A = 0.80*KCP_N + 0.2*KCP_O
RT = 298.0 * 0.0821 # atm.L.mol-1

# find_x fnction return water length in a vessel by using Lv, Tc and Rdv
def find_x(l,ct,bs):
    x = 0
    count = 0
    if (bs ** 2 / 0.127 ** 2 -1) * 1000 * ct + 100 + cp > BP:
        while 1:
            count = count + 1
            if(count >=20):
                x = l
                break
            fr = x / l
            fwt = 1.0 - center_PLC/100.0/10.0
            fbt = center_PLC/100.0/10.0
            fwe = fr*fbt + fwt
            fbe = fbt*(1.0-fr)
            funcKt = KCP_A*fwt + fbt/RT
            funcKe = KCP_A*fwe + fbe/RT
            EquilP = BP * funcKt/funcKe
            # judge here is Pw+Pc-Pb
            judge = ((abs(bs) - x) ** 2 / 0.127 ** 2 -1) * ct * 1000.0 + 100.0 + cp - EquilP
            # slope here is the slope of "judge" at x
            slope = - ct * 1000.0 / 0.127 ** 2 * 2.0 * (abs(bs) -x) + BP*funcKt/(funcKe**2) * (fbt/l*KCP_A - fbt/RT/l)
            # Newton Iteration
            x = x - judge / slope
            #print judge,slope,x
            if (abs(judge) < 0.0001):
                # To judge when to stop
                break
            # A Statement to avoid overflow
            if x >= l or x < 0:
                x = l/1.2
    return x

# A function to compute the resistance of the embolized vessels at any given slice
def drawf_PLC(ct,vl):
    # Begin the slice array from the distal end of region a by the thickness of 0.001 m
    site = -0.137
    dsite = 0.001
    # Stop after the last slice at 0.137 m
    while(site < 0.138):
        # The resistance of slices in a/e is constant R0, we use R0=1.0
        if(abs(site) > 0.127):
            NEWR.append(1.0)
            OLDR.append(1.0)
        # The resistance of slices of end and normal type
        elif(abs(site) >= vl):
            klist = []
            ith = 1
            while(ith <= 100):
                tmpml = ith*0.01*vl + 0.127-abs(site)
                if(tmpml > vl):
                    tmpml = vl
                    tmpbs = abs(site) - (vl - ith*0.01*vl)
                    tmpwl = find_x(tmpml,ct,tmpbs)
                    tmpki = 2.0*tmpml*tmpwl / (tmpml**2 + tmpwl**2)
                    klist.append(tmpki)
                else:
                    tmpwl = find_x(tmpml,ct,0.127)
                    tmpki = 2.0*tmpml*tmpwl / (tmpml**2 + tmpwl**2)
                    klist.append(tmpki)
                ith += 1
            kr = numpy.mean(klist)
            NRdx = 100.0 / ((100.0-center_PLC) + kr*center_PLC)
            ORdx = 100.0 / (100.0-center_PLC)
            NEWR.append(NRdx)
            OLDR.append(ORdx)
        # The resistance of slices in region c
        else:
            klist = []
            ith = 1
            while(ith <= 100):
                tmpml = vl
                tmpbs = max(abs(site)+vl-ith*0.01*vl,abs(abs(site)-ith*0.01*vl))
                tmpes = min(abs(site)+vl-ith*0.01*vl,abs(abs(site)-ith*0.01*vl))
                if(abs(tmpbs) < vl):
                    # A judgement use to judge whether the equilibrium establish from both ends
                    ljudge = vl / (2*tmpes) * BP
                    # Patm = 100.0, Tension in MPa
                    rjudge = cp + ((tmpes/0.127)**2 - 1.0)*ct*1000.0 + 100.0
                    if(ljudge >= rjudge):
                        # When the slice is normal type slice
                        tmpwl = find_x(vl,ct,tmpbs)
                        tmpki = 2.0*tmpml*tmpwl / (tmpml**2 + tmpwl**2)
                        klist.append(tmpki)
                    else:
                        # When the slice is center type slice
                        tmpwl = 2.0 * find_x(vl/2.0,ct,vl/2.0)
                        tmpki = 2.0*tmpml*tmpwl / (tmpml**2 + tmpwl**2)
                        klist.append(tmpki)
                else:
                    tmpwl = find_x(vl,ct,tmpbs)
                    tmpki = 2.0*tmpml*tmpwl / (tmpml**2 + tmpwl**2)
                    klist.append(tmpki)
                ith += 1
            kr = numpy.mean(klist)
            NRdx = 100.0 / ((100.0-center_PLC) + kr*center_PLC)
            ORdx = 100.0 / (100.0-center_PLC)
            NEWR.append(NRdx)
            OLDR.append(ORdx)
        site += dsite

# main part
kmaxs = [2.38,6.33,5.02,4.28,2.35,5.95,1.29,3.28,4.69,1.57,2.74,5.37,5.36,
         9.61,8.95,13.7,4.56,12.1,7.67,5.78,12.1,5.68,15.6,7.07,8.40,6.74,
         8.07,9.86,11.4,7.88,5.51,9.13,22.3,10.8,5.25]
cpres = [7.8023]*13 + [6.2637]*22
vlens = [0.02877]*13 + [0.07105]*22
ktens = [1.28,3.14,2.25,1.64,.829,1.22,.876,1.04,1.93,.956,2.35,4.34,4.24,
         4.92,5.78,1.36,3.32,2.67,3.81,4.52,5.86,3.73,4.04,3.75,1.84,2.75,
         4.96,2.35,2.60,3.72,2.71,4.09,5.78,4.86,2.78]
k800s = [1.55,3.37,2.90,1.85,1.13,1.84,1.15,1.55,2.15,1.34,2.59,4.86,4.82,
         6.58,8.39,4.05,4.42,4.81,5.48,5.59,7.07,4.59,5.58,4.16,2.11,3.45,
         7.08,3.14,3.46,5.24,3.69,5.64,7.09,8.21,3.12]
cp = 6.2637
vessel_length = 0.07105
#cp = 7.8023
#vessel_length = 0.02877
for n in range(len(kmaxs)):
    tmpkmax = kmaxs[n]
    tmpkten = ktens[n]
    tmpk800 = k800s[n]
    center_PLC = (100.0 - 100.0*tmpkten/tmpkmax) * 274.0/254.0
    maxBP = 100.0
    minBP = 0.0
    BP = (maxBP+minBP)/2.0
    while(1):
        OLDR = []
        NEWR = []
        drawf_PLC(0.0565,vessel_length)
        tempkvalue = 275.0/sum(NEWR) * tmpkmax
        print tempkvalue,tmpk800,BP
        if(abs(tempkvalue/tmpk800 - 1.0) > 0.0001):
            if(tempkvalue > tmpk800):
                minBP = BP
                BP = (maxBP+minBP)/2.0
                print BP
            else:
                maxBP = BP
                BP = (maxBP+minBP)/2.0
                print BP
        else:
            break
    judge = raw_input("Continue (Y?N) >>")
    if(judge == "N" or judge == "n"):
        break
    
    
