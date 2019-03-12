#Run FRG
import config as cfg 
import numpy as np
npatch = cfg.npatch
qsc = 0
qsdw = cfg.pi 
qcdw = cfg.pi 
V = np.zeros((npatch,npatch))
Vsc,Vsdw,Vcdw = 0,0,0


def lindhard(La, Ek1, vK1, Ek2, vK2, radius):
	z1 = (1j * La - Ek1)/ vK1
	z2 = (1j * La - Ek2)/ vK2
	if abs(z1-z2)<1.0e-10 :
		return 2./cfg.twopi * np.real(-2 * radius / (radius * radius - z1 * z1) / ( vK1 * vK1 ))
	else:
		return 2./cfg.twopi * np.real((np.log( (radius-z1)/(-radius-z1) ) - np.log( (radius-z2)/(-radius-z2) ) ) / (z1-z2) / ( vK1 * vK2 ))
def runfrg(Lc,kc):

	#>>>>>>>>>>>>>>>>---calculate bubbles--->>>>>>>>>>>>>>>>>
	for i5 in range(npatch):
		for i6 in range(npatch):
			cfg.xpp[i5,i6] = -lindhard(Lc,cfg.FermiE[i5], cfg.FermiV[i5],\
				-cfg.FermiE[i6], cfg.FermiV[i6],kc)
			cfg.xph[i5,i6] = lindhard(Lc,cfg.FermiE[i5], cfg.FermiV[i5],\
				cfg.FermiE[i6], cfg.FermiV[i6],kc)
	#<<<<<<<<<<<<<<<<---calculate bubbles----<<<<<<<<<<<<<<<<<<


	cfg.dGa = np.zeros((npatch,npatch,npatch))
	for k3 in range(npatch):
		for k2 in range(npatch):
			for k1 in range(npatch):
				k4 = cfg.leg4[k1,k2,k3]
				if k4 < 0 :
					continue
				
				for k5 in range(npatch):
					#Particle - Hole channel
					k6 = cfg.leg4[k1,k2,k5]
					if k6 != -1:
						cfg.dGa[k1, k2, k3] = cfg.dGa[k1, k2, k3] \
						+ cfg.Ga[k1, k2, k5] * cfg.xpp[k5, k6] * cfg.Ga[k6, k5, k3]
					#Particle - Particle cross channel
					k6 = cfg.leg4[k1, k5, k3]
					if k6 != -1:
						cfg.dGa[k1, k2, k3] = cfg.dGa[k1, k2, k3] \
						+ cfg.Ga[k1, k5, k3] * cfg.xph[k5, k6] * cfg.Ga[k6, k2, k5]
					
					#Particle - Particle direct channel
					k6 = cfg.leg4[k5, k2, k3]
					if k6 != -1:
						cfg.dGa[k1, k2, k3] = cfg.dGa[k1, k2, k3] \
						- 2 * cfg.Ga[k1, k6, k5] * cfg.xph[k5, k6] * cfg.Ga[k5, k2, k3]

					k6 = cfg.leg4[k5, k2, k3]
					if k6 != -1:
						cfg.dGa[k1, k2, k3] = cfg.dGa[k1, k2, k3] \
						+ cfg.Ga[k1, k6, k4] * cfg.xph[k5, k6] * cfg.Ga[k5, k2, k3]

					k6 = cfg.leg4[k2, k5, k3]
					if k6 != -1:
						cfg.dGa[k1, k2, k3] = cfg.dGa[k1, k2, k3] \
						+ cfg.Ga[k1, k6, k5] * cfg.xph[k5, k6] * cfg.Ga[k2, k5, k3]

#------------small functions----------
def findpatch(i,Q,info):
	for k in range(npatch):
		if abs((cfg.FermiM[i] + info * cfg.FermiM[k] - Q) % cfg.twopi) <1.0e-5 :
			return k
		else:
			return -1
def check(ndim,V):
	for i in range(ndim):
		for j in range(i-1):
			if abs(V[i,j] - V[j,i]) > 1.e-6 :
				print('hermiticity violated')
#-------------------------------------

def pdw():
	V = np.zeros((npatch,npatch))
	for i1 in range(npatch):
		i2 = findpatch(i1, qsc, 1)
		if i2 < 0 : continue
		for i4 in range(npatch):
			i3 = findpatch(i4, qsc, 1)
			if i3 < 0: continue
			V[i1,i4] = V[i1,i4] + cfg.Ga[i1,i2,i3]
	V = V/npatch
	check(npatch,V)
	V = -V
	return sorted(abs(np.linalg.eigvals(V)))[-1]

def sdw():
	V = np.zeros((npatch,npatch))
	for i1 in range(npatch):
		i3 = findpatch(i1, qsdw, -1)
		if i3 < 0 : continue
		for i4 in range(npatch):
			i2 = findpatch(i4, qsdw, -1)
			if i2 < 0: continue
			V[i1,i4] = V[i1,i4] + cfg.Ga[i1,i2,i3]
	V = V/npatch
	check(npatch,V)
	return sorted(abs(np.linalg.eigvals(V)))[-1]

def cdw():
	V = np.zeros((npatch,npatch))
	for i1 in range(npatch):
		i4 = findpatch(i1, qcdw, -1)
		if i4 < 0 : continue
		for i3 in range(npatch):
			i2 = findpatch(i3, qcdw, -1)
			if i2 < 0: continue
			V[i1,i3] = V[i1,i3] - cfg.Ga[i1,i2,i3] * 2

	for i1 in range(npatch):
		i3 = findpatch(i1, qcdw, -1)
		if i3 < 0 : continue
		for i4 in range(npatch):
			i2 = findpatch(i4, qcdw, -1)
			if i2 < 0: continue
			V[i1,i4] = V[i1,i4] + cfg.Ga[i1,i2,i3]

	V = V/npatch
	check(npatch,V)
	return sorted(abs(np.linalg.eigvals(V)))[-1]


def analysis(Lc,scheme):
	if scheme == "current" :
		mi, mj, mk = 1+2*cfg.extra, 2+4*cfg.extra, 3+6*cfg.extra

		b11r=2*cfg.Ga[mj,mk,mk]-4*cfg.Ga[mk,mj,mk]
		b11s=-cfg.Ga[mj,mk,mk]-1/2*cfg.Ga[mk,mj,mk]+1/2*(-cfg.Ga[mj,mk,mk]+cfg.Ga[mk,mj,mk])
		b12r=2*cfg.Ga[0,mi,mk]-4*cfg.Ga[mi,0,mk]
		b12s=-cfg.Ga[0,mi,mk]-1/2*cfg.Ga[mi,0,mk]+1/2*(-cfg.Ga[0,mi,mk]+cfg.Ga[mi,0,mk])
		f12r=2*cfg.Ga[0,mk,mk]-4*cfg.Ga[mk,0,mk]
		f12s=-cfg.Ga[0,mk,mk]-1/2*cfg.Ga[mk,0,mk]+1/2*(-cfg.Ga[0,mk,mk]+cfg.Ga[mk,0,mk])
		u11r=4*cfg.Ga[0,0,mk]+4*cfg.Ga[mk,mk,0]
		u12r=2*cfg.Ga[0,mj,mi]+2*cfg.Ga[mi,mk,mj]+2*cfg.Ga[mj,0,mi]+2*cfg.Ga[mk,mi,mj]
		u12s=-cfg.Ga[0,mj,mi]-cfg.Ga[mk,mi,mj]
		
		data = str([Lc, 1/abs(b11r), 1/abs(b11s), 1/abs(b12r), \
		1/abs(b12s), 1/abs(f12r), 1/abs(f12s), 1/abs(u11r), 1/abs(u12r), 1/abs(u12s)])+"\n"
		data = data.replace("[","")
		data = data.replace("]","")
		return data
	elif scheme == "channel":
		Vsc = pdw()
		Vsdw = sdw()
		Vcdw = cdw()
		data = str([Lc, Vsc, Vsdw, Vcdw, 1/abs(Vsc), 1/abs(Vsdw), 1/abs(Vcdw)])+"\n"
		data = data.replace("[","")
		data = data.replace("]","")
		return data


