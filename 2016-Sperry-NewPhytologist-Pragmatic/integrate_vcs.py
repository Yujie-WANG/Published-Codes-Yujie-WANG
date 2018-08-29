from numpy import array
import os

target_p = []
target_e = []

rfile = open("./output/multi_vc_e.txt", "r")
for i in rfile:
	tmp = map(float, i.split("\t"))
	target_p.append(tmp[0])
	target_e.append(tmp[1])
rfile.close()

k_target = sum(target_e)
m_target = target_e[0]

fb = 2.0
fc = 5.0
judge = 0
while(1):
	print "\n",fb,fc,"\n"
	cmd_str = "integrate_vcs_fit " + str(fb) + " " + str(fc) + " 1.3 3.0"
	os.system(cmd_str)
	fit_p = []
	fit_e = []
	tmp_file = open("./output/single_vc_e.txt","r")
	for i in tmp_file:
		tmp = map(float, i.split("\t"))
		fit_p.append(tmp[0])
		fit_e.append(tmp[1])
	tmp_file.close()
	k_fit = sum(fit_e)
	m_fit = fit_e[0]
	print "\nTarget-fit are:",m_target-m_fit,k_target-k_fit
	if(abs(m_target - m_fit) < 1E-4 and abs(k_target - k_fit) < 5E-2):
		break
	if(judge == 0):
		fbb = m_target / m_fit * fb
		fb = fbb
		judge = 1
	else:
		fcc = k_fit / k_target * fc
		fc = fcc
		judge = 0
