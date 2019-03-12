#main program
import config as cfg 
import bandinfo_fun as bandinfo 
import runfrg as runfrg
import smfrg_new as smfrg_new
import numpy as np
import time

bandinfo.setbandinfo()


scheme = "SMFRG_new"

print(cfg.qrange,'\n',cfg.R)

wIR, wUV, divergence = 1.0e-6, 100, 60
La = wUV
b = np.exp(np.log(wUV/wIR) / 600)
rc = bandinfo.kspan/bandinfo.np4


while La > wIR :
#for t in range(1):
#	print(La, end="\n")
	print( time.asctime(time.localtime(time.time())))
	Lc, dL = La * (1 + 1/b)/2, La - La/b
	La = La / b

	if scheme == "FRG" :
		runfrg.runfrg(Lc,rc)
		cfg.dGa = np.array(cfg.dGa[:,:,:]) * dL / cfg.twopi
		cfg.Ga = np.add(cfg.Ga, cfg.dGa)
		#if cfg.Ga.max() > divergence:
		#	break
		if runfrg.Vsc > divergence or runfrg.Vsdw > divergence or runfrg.Vcdw > divergence:
			break
		data = runfrg.analysis(Lc,"channel")
		print(data)
		with open('flowdata.txt', 'a') as f:
			f.write(data)



	elif scheme == "SMFRG_new":
		smfrg_new.smfrg(Lc,dL)


		data = smfrg_new.analysis_sm(Lc)
		print(data[0])
		if abs(data[1]) > divergence or abs(data[2]) > divergence or abs(data[3]) > divergence:
			break
		with open('flowdata.txt', 'a') as f:
			f.write(data[0])
		
	else:
		break


