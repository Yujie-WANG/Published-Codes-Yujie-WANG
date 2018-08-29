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
def residual(p,plsqt,plsqk,plsqkmax):
    global BP, center_PLC
    plsqbp,plsqplc = p
    BP = plsqbp
    center_PLC = plsqplc
    #print BP, center_PLC
    tmpk = []
    for tension in plsqt:
        print f,tension
        global NEWR, OLDR
        NEWR = []
        OLDR = []
        drawf_PLC(tension,vessel_length)
        tempvalue = 275.0/sum(NEWR) * plsqkmax
        tmpk.append(tempvalue)
    results = plsqk - tmpk
    print results/plsqk * 100
    return results


# center_tension and k_n
p1 = open("./Populus-Curve-Data/Populus.PLC.BP.1.result.txt")
p2 = open("./Populus-Curve-Data/Populus.PLC.BP.2.result.txt")
p3 = open("./Populus-Curve-Data/Populus.PLC.BP.3.result.txt")
p4 = open("./Populus-Curve-Data/Populus.PLC.BP.4.result.txt")
p5 = open("./Populus-Curve-Data/Populus.PLC.BP.5.result.txt")
p6 = open("./Populus-Curve-Data/Populus.PLC.BP.6.result.txt")
p7 = open("./Populus-Curve-Data/Populus.PLC.BP.7.result.txt")
p8 = open("./Populus-Curve-Data/Populus.PLC.BP.8.result.txt")
p9 = open("./Populus-Curve-Data/Populus.PLC.BP.9.result.txt")
p10 = open("./Populus-Curve-Data/Populus.PLC.BP.10.result.txt")
p11 = open("./Populus-Curve-Data/Populus.PLC.BP.11.result.txt")
p12 = open("./Populus-Curve-Data/Populus.PLC.BP.12.result.txt")
p13 = open("./Populus-Curve-Data/Populus.PLC.BP.13.result.txt")
p14 = open("./Populus-Curve-Data/Populus.PLC.BP.14.result.txt")
p15 = open("./Populus-Curve-Data/Populus.PLC.BP.15.result.txt")
p16 = open("./Populus-Curve-Data/Populus.PLC.BP.16.result.txt")
p17 = open("./Populus-Curve-Data/Populus.PLC.BP.17.result.txt")
p18 = open("./Populus-Curve-Data/Populus.PLC.BP.18.result.txt")
p19 = open("./Populus-Curve-Data/Populus.PLC.BP.19.result.txt")
p20 = open("./Populus-Curve-Data/Populus.PLC.BP.20.result.txt")
p21 = open("./Populus-Curve-Data/Populus.PLC.BP.21.result.txt")
p22 = open("./Populus-Curve-Data/Populus.PLC.BP.22.result.txt")
a4 = open("./Acer-Curve-Data/Acer.PLC.BP.4.result.txt")
a5 = open("./Acer-Curve-Data/Acer.PLC.BP.5.result.txt")
a6 = open("./Acer-Curve-Data/Acer.PLC.BP.6.result.txt")
a7 = open("./Acer-Curve-Data/Acer.PLC.BP.7.result.txt")
a8 = open("./Acer-Curve-Data/Acer.PLC.BP.8.result.txt")
a9 = open("./Acer-Curve-Data/Acer.PLC.BP.9.result.txt")
a10 = open("./Acer-Curve-Data/Acer.PLC.BP.10.result.txt")
filenames = (p1,p2,p3,p4,p5,p6,p7,p8,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,
             a4,a5,a6,a7,
             a8,a9,a10)
kmaxs = [7.87792820555e-05,5.51111059804e-05,9.1308307757e-05,0.000223140158687,5.2513798422e-05,
         0.000107552209141,9.61041527914e-05,8.94575573947e-05,4.56364869361e-05,
         0.000121091002867,7.66738220636e-05,5.78E-05,0.000121196260621,5.67852519865e-05,
         0.000155794096671,7.07272566289e-05,8.40455646959e-05,6.73889527553e-05,8.06517725917e-05,
         9.86288576902e-05,0.000113691214531,
         5.02241952875e-05,4.28454330251e-05,5.96340396812e-05,1.29161075758e-05,
         3.27923081379e-05,4.68507521228e-05,1.56919028048e-05]
cps = [6.2637,6.2637,6.2637,6.2637,6.2637,6.2637,6.2637,6.2637,6.2637,
       6.2637,6.2637,6.2637,6.2637,6.2637,6.2637,6.2637,6.2637,6.2637,6.2637,
       6.2637,6.2637,
       7.8023,7.8023,7.8023,7.8023,
       7.8023,7.8023,7.8023]
vls = [0.07105,0.07105,0.07105,0.07105,0.07105,0.07105,0.07105,0.07105,0.07105,
       0.07105,0.07105,0.07105,0.07105,0.07105,0.07105,0.07105,0.07105,0.07105,0.07105,
       0.07105,0.07105,
       0.02877,0.02877,0.02877,0.02877,
       0.02877,0.02877,0.02877]


savefile = open("./Fitting-Results.txt","a")
bplcfile = open("./BP-PLC.txt","a")
for i in range(len(filenames)):
    f = filenames[i]
    cp = cps[i]
    vessel_length = vls[i]
    if(f):
        fileresults = []
        tensions = []
        kh25s = []
        tmpkmax = kmaxs[i]
        print "Success in openning file."
        for j in f:
            fileresults.append(map(float,j.split('\t')))
        for g in fileresults:
            tensions.append(g[0])
            kh25s.append(g[3])
        p0 = [60,50]
        Sigt = numpy.array(tensions)
        Sigk = numpy.array(kh25s)
        print Sigt
        print Sigk
        plsq = leastsq(residual,p0,args=(Sigt,Sigk,tmpkmax))
        print "The best fitting of [BP,PLC] and Kmax is"
        print plsq[0],tmpkmax
        plsq_results = residual(plsq[0],Sigt,Sigk,tmpkmax)
        print plsq_results
        tempsum = 0.0
        for e in range(len(plsq_results)):
            tempsum += abs(plsq_results[e]/kh25s[e])
        tempave = tempsum / (len(plsq_results)+1.0)
        print "\n\n\n"
        # Kmax, VL, BP, PLC
        bplcfile.write(str(tmpkmax)+"\t"+str(vessel_length)+"\t"+str(plsq[0][0])+"\t"+str(plsq[0][1])+"\t"+str(tempave)+"\n")
        savefile.write(str(f)+"\n")
        savefile.write("The best fitting of [BP,PLC] and Kmax is:\n")
        savefile.write(str(plsq[0])+str(tmpkmax)+"\n")
        savefile.write("Tension\tKh25s\tK-Kpre\n")
        for i in range(len(Sigk)):
            savefile.write(str(Sigt[i])+"\t"+str(Sigk[i])+"\t"+str(plsq_results[i])+"\n")
        savefile.write("\n\n")
        f.close()
    else:
        print "Wrong File name!"
savefile.close()
bplcfile.close()
