import numpy

KCP_O = 1.30E-3 #mol.L-3.atm-1
KCP_N = 6.10E-4
KCP_A = 0.80*KCP_N + 0.2*KCP_O
RT = 298.0 * 0.0821 # atm.L.mol-1

# range of central tension from 0.0 to 5.0 MPa
center_tension = []
tmptension = 0.0
while(tmptension <= 5.0):
    center_tension.append(tmptension)
    tmptension += 0.01

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
            
# Main Body of the Code
# cp here is actually capillary pressure
# Here we use lists to store the values of each parameter

cps = [12.0,12.0,12.0,12.0,12.0]
vessel_lengths = [0.05,0.05,0.05,0.05,0.05]
filenames = ["./BP-curve/BP-10.txt","./BP-curve/BP-30.txt","./BP-curve/BP-50.txt","./BP-curve/BP-70.txt","./BP-curve/BP-90.txt"]
Kmaxs = [1.0E-4,1.0E-4,1.0E-4,1.0E-4,1.0E-4]
cPLCs = [50.0,50.0,50.0,50.0,50.0]
BPs = [10.0,30.0,50.0,70.0,90.0]
# main part
for i in range(len(filenames)):
    # Assign the parameters!
    OLDR = []
    NEWR = []
    savefilename = filenames[i]
    kmax = Kmaxs[i]
    center_PLC = cPLCs[i]
    BP = BPs[i]
    cp = cps[i]
    vessel_length = vessel_lengths[i]
    savefile = open(savefilename,"w+")
    savefile.write("Tension\tKh\tPLC\n")
    # Here we give a limit to vessel length to avoid overflow
    if vessel_length < 0.127:
    #cycle through center tensions
        for tension in center_tension:
        # Clear OLDR,NEWR in each cycle
            OLDR = []
            NEWR = []
            # Compute the resistance in 275 slices
            drawf_PLC(tension,vessel_length)
            tmpkh = kmax * (274.0 / sum(NEWR))
            tmplc = 100.0 - 100.0*(274 / sum(NEWR))
            tmpstr = str(tension) + "\t" + str(tmpkh) + "\t" + str(tmplc)
            print tmpstr
            savefile.write(tmpstr + "\n")
    # Save the output results and close the file.
    savefile.close()