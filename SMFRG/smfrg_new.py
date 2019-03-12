#smfrg
import config as cfg 
import numpy as np 
import bandinfo_fun as bi
import scipy
from scipy.integrate import quad
from bandinfo_fun import Gf
from config import fm,fmc
from config import nq, nfactor, nsublattice, NNN, nR
from config import dPchanR, dCchanR, dDchanR
from config import fourierscheme, bubblescheme


nf = nfactor
ns = nsublattice
dPchanRs, dCchanRs,dDchanRs = np.zeros((nR,ns,ns,nf,ns,ns,nf),dtype=complex),\
np.zeros((nR,ns,ns,nf,ns,ns,nf),dtype=complex), np.zeros((nR,ns,ns,nf,ns,ns,nf),dtype=complex)


#fm(k)[a][b][m]  form factor, a,b are sublattice index
#other label by [nR/nq,(abm),(cdn)]

################################################define bubble integration


def xpp_fun(La, p, q, m, n, a, b, c, d):
	res = fmc(p)[a][b][m] * Gf(p+q,+1,La)[a][d]\
	* Gf(-p,-1,La)[b][c] * fm(p)[c][d][n]
	return np.real(-res/cfg.pi)

def xph_fun(La, p, q, m, n, a, b, c, d):
	res = fmc(p)[a][b][m] * Gf(p+q,+1,La)[d][a]\
	* Gf(p,+1,La)[b][c] * fm(p)[c][d][n]
	return np.real(-res/cfg.pi)


def xpp_int(La, q, m, n, a, b, c, d, flag):
	if flag == 'integration':
		if La > 0.01:
			result, err = quad(lambda p:\
			xpp_fun(La, p, q, m, n, a, b, c, d), -cfg.pi,cfg.pi)
			return [np.real(result), err]
		else:
			dkk = 0.1
			kk0=[-np.pi, -np.arccos(bi.mu/2) - dkk, -np.arccos(bi.mu/2) + dkk, np.arccos(bi.mu/2) - dkk, np.arccos(bi.mu/2) + dkk, \
			-q-np.arccos(bi.mu/2) - dkk, -q-np.arccos(bi.mu/2) + dkk, -q+np.arccos(bi.mu/2) - dkk, -q+np.arccos(bi.mu/2) + dkk, np.pi]
			for kk in range(len(kk0)):
				if kk0[kk] > np.pi:
					kk0[kk] -= 2*np.pi
				elif kk0[kk] < -np.pi:
					kk0[kk] += 2*np.pi
			kk0 = sorted(set(kk0))
			res = 0
			summ = 0
			summe = 0
			for ii in range(len(kk0)-1):
				res, err = quad(lambda p:\
			xpp_fun(La, p, q, m, n, a, b, c, d), kk0[ii], kk0[ii+1])
				summ += res
				summe += err
			return [np.real(summ), summe]

	elif flag == 'sum':
		res = 0
		for pi in range(nq-1):
			res += xpp_fun(La, cfg.qrange[pi], q, m, n, a, b, c, d) / (nq-1) * cfg.twopi
		return np.real(res)



def xph_int(La, q, m, n, a, b, c, d, flag):
	if flag == 'integration':
		if La > 0.01:
			result, err = quad(lambda p:\
			xph_fun(La, p, q, m, n, a, b, c, d), -cfg.pi,cfg.pi)
			return [np.real(result), err]
		else:
			dkk = 0.1
			kk0=[-np.pi, -np.arccos(bi.mu/2) - dkk, -np.arccos(bi.mu/2) + dkk, np.arccos(bi.mu/2) - dkk, np.arccos(bi.mu/2) + dkk, \
			-q-np.arccos(bi.mu/2) - dkk, -q-np.arccos(bi.mu/2) + dkk, -q+np.arccos(bi.mu/2) - dkk, -q+np.arccos(bi.mu/2) + dkk, np.pi]
			for kk in range(len(kk0)):
				if kk0[kk] > np.pi:
					kk0[kk] -= 2*np.pi
				elif kk0[kk] < -np.pi:
					kk0[kk] += 2*np.pi
			kk0 = sorted(set(kk0))
			res = 0
			summ = 0
			summe = 0
			for ii in range(len(kk0)-1):
				res, err = quad(lambda p:\
			xph_fun(La, p, q, m, n, a, b, c, d), kk0[ii], kk0[ii+1])
				summe += err
				summ += res
			return [np.real(summ), summe]
	elif flag == 'sum':
		res = 0
		for p in range(nq-1):
			res += xph_fun(La, cfg.qrange[p], q, m, n, a, b, c, d) / (nq-1) * cfg.twopi
		return np.real(res)
	




def FourierTrans():
	cfg.dPchanR, cfg.dCchanR, cfg.dDchanR = np.zeros((cfg.nR,NNN,NNN),dtype=complex),\
	np.zeros((cfg.nR,NNN,NNN),dtype=complex), np.zeros((cfg.nR,NNN,NNN),dtype=complex)

	for m in range(NNN): 
		for n in range(NNN):
			for Ri in range(cfg.nR): # tabulate m,n,R
			#--------------------------------------
				for q in range(cfg.nq-1): #sum over q
					cfg.dPchanR[Ri][m][n] = cfg.dPchanR[Ri][m][n] +\
					cfg.dPchan[q][m][n] * np.exp(1j*cfg.qrange[q]*cfg.R[Ri])/ (cfg.nq - 1)

					cfg.dCchanR[Ri][m][n] = cfg.dCchanR[Ri][m][n] +\
					cfg.dCchan[q][m][n] * np.exp(1j*cfg.qrange[q]*cfg.R[Ri])/ (cfg.nq - 1)

					cfg.dDchanR[Ri][m][n] = cfg.dDchanR[Ri][m][n] +\
					cfg.dDchan[q][m][n] * np.exp(1j*cfg.qrange[q]*cfg.R[Ri])/ (cfg.nq - 1)
	

def InvFourierTrans():
	cfg.dPchan, cfg.dCchan, cfg.dDchan = np.zeros((cfg.nq,NNN,NNN)),\
	np.zeros((cfg.nq,NNN,NNN)), np.zeros((cfg.nq,NNN,NNN))

	for m in range(NNN): 
		for n in range(NNN):
			for q in range(nq): # tabulate m,n,q
			#--------------------------------------
				for Ri in range(nR-1): #sum over R
					cfg.dPchan[q][m][n] = np.real(cfg.dPchan[q][m][n] +\
					cfg.dPchanRQ[Ri][m][n] * np.exp(-1j*cfg.qrange[q]*cfg.R[Ri]))

					cfg.dCchan[q][m][n] = np.real(cfg.dCchan[q][m][n] +\
					cfg.dCchanRQ[Ri][m][n] * np.exp(-1j*cfg.qrange[q]*cfg.R[Ri]))
	
					cfg.dDchan[q][m][n] = np.real(cfg.dDchan[q][m][n] +\
					cfg.dDchanRQ[Ri][m][n] * np.exp(-1j*cfg.qrange[q]*cfg.R[Ri]))


def smfrg(Lc, dL):
	nR2 = int((nR-1)/2)
	cfg.xpps = np.reshape(cfg.xpps,(nq,ns,ns,nf,ns,ns,nf))
	cfg.xphs = np.reshape(cfg.xphs,(nq,ns,ns,nf,ns,ns,nf))
	
	for q in range(nq):
		for a in range(ns):
			for b in range(ns):
				for m in range(nf):
					for c in range(ns):
						for d in range(ns):
							for n in range(nf):
								cfg.xpps[q,a,b,m,c,d,n] = \
								xpp_int(Lc, cfg.qrange[q], m, n,a,b,c,d,'integration')[0]

								cfg.xphs[q,a,b,m,c,d,n] = \
								xph_int(Lc, cfg.qrange[q], m, n,a,b,c,d,'integration')[0]

	cfg.xpps = np.reshape(cfg.xpps,(nq,NNN,NNN))
	cfg.xphs = np.reshape(cfg.xphs,(nq,NNN,NNN))
	#increasement of each channel
	cfg.dPchan, cfg.dCchan, cfg.dDchan = np.zeros((cfg.nq,NNN,NNN)),\
	np.zeros((cfg.nq,NNN,NNN)), np.zeros((cfg.nq,NNN,NNN))

	for q in range(nq):	
		cfg.dPchan[q] = np.dot(np.dot(cfg.Pchan[q],cfg.xpps[q]),cfg.Pchan[q])
		cfg.dCchan[q] = np.dot(np.dot(cfg.Cchan[q],cfg.xphs[q]),cfg.Cchan[q])
		cfg.dDchan[q] = (np.dot(np.dot((cfg.Cchan[q] - cfg.Dchan[q]),cfg.xphs[q]),cfg.Dchan[q])\
		+ np.dot(np.dot(cfg.Dchan[q],cfg.xphs[q]),(cfg.Cchan[q] - cfg.Dchan[q])))


	#---------------------Contact each channel---------------
	cfg.dPchanR, cfg.dCchanR, cfg.dDchanR = np.zeros((cfg.nR,NNN,NNN),dtype=complex),\
	np.zeros((cfg.nR,NNN,NNN),dtype=complex), np.zeros((cfg.nR,NNN,NNN),dtype=complex)

	dPchanRs, dCchanRs, dDchanRs = np.zeros((cfg.nR,ns,ns,nf,ns,ns,nf),dtype=complex),\
	np.zeros((cfg.nR,ns,ns,nf,ns,ns,nf),dtype=complex), np.zeros((cfg.nR,ns,ns,nf,ns,ns,nf),dtype=complex)
	cfg.dPchanRQ, cfg.dCchanRQ, cfg.dDchanRQ = np.zeros((cfg.nR,ns,ns,nf,ns,ns,nf),dtype=complex),\
	np.zeros((cfg.nR,ns,ns,nf,ns,ns,nf),dtype=complex), np.zeros((cfg.nR,ns,ns,nf,ns,ns,nf),dtype=complex)

	FourierTrans()
	

	##in real space
	dPchanRs = np.reshape(cfg.dPchanR,(nR,ns,ns,nf,ns,ns,nf))
	dCchanRs = np.reshape(cfg.dCchanR,(nR,ns,ns,nf,ns,ns,nf))
	dDchanRs = np.reshape(cfg.dDchanR,(nR,ns,ns,nf,ns,ns,nf))

	#################

	#---------------------------check------------------

	for ii in range(cfg.nR):
		dPchanRs[ii] = [[[[[[dPchanRs[ii, 0, 0, 0, 0, 0, 0], dPchanRs[ii, 0, 0, 0, 0, 0, 1], dPchanRs[ii, 0, 0, 0, 0, 0, 2]]]], [[[dPchanRs[ii, 0, 0, 1, 0, 0, 0], dPchanRs[ii, 0, 0, 1, 0, 0, 1], dPchanRs[ii, 0, 0, 1, 0, 0, 2]]]], [[[dPchanRs[ii, 0, 0, 2, 0, 0, 0], dPchanRs[ii, 0, 0, 2, 0, 0, 1], dPchanRs[ii, 0, 0, 2, 0, 0, 2]]]]]]]
		dCchanRs[ii] = [[[[[[dCchanRs[ii, 0, 0, 0, 0, 0, 0], dCchanRs[ii, 0, 0, 0, 0, 0, 1], dCchanRs[ii, 0, 0, 0, 0, 0, 2]]]], [[[dCchanRs[ii, 0, 0, 1, 0, 0, 0], dCchanRs[ii, 0, 0, 1, 0, 0, 1], dCchanRs[ii, 0, 0, 1, 0, 0, 2]]]], [[[dCchanRs[ii, 0, 0, 2, 0, 0, 0], dCchanRs[ii, 0, 0, 2, 0, 0, 1], dCchanRs[ii, 0, 0, 2, 0, 0, 2]]]]]]]
		dDchanRs[ii] = [[[[[[dDchanRs[ii, 0, 0, 0, 0, 0, 0], dDchanRs[ii, 0, 0, 0, 0, 0, 1], dDchanRs[ii, 0, 0, 0, 0, 0, 2]]]], [[[dDchanRs[ii, 0, 0, 1, 0, 0, 0], dDchanRs[ii, 0, 0, 1, 0, 0, 1], dDchanRs[ii, 0, 0, 1, 0, 0, 2]]]], [[[dDchanRs[ii, 0, 0, 2, 0, 0, 0], dDchanRs[ii, 0, 0, 2, 0, 0, 1], dDchanRs[ii, 0, 0, 2, 0, 0, 2]]]]]]]

	#---------------------------check------------------
		

	cfg.dPchanRQ, cfg.dCchanRQ, cfg.dDchanRQ = dPchanRs, dCchanRs, dDchanRs
	print(cfg.dPchanRQ)
	#---------------------------------projection--------------------------
	nR2 = int((nR-1)/2)
	R = -2
	cfg.dPchanRQ[R + nR2] = [[[[[[dPchanRs[-2 + nR2, 0, 0, 0, 0, 0, 0], dPchanRs[-2 + nR2, 0, 0, 0, 0, 0, 1], dPchanRs[-2 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dPchanRs[-2 + nR2, 0, 0, 1, 0, 0, 0], dPchanRs[-2 + nR2, 0, 0, 1, 0, 0, 1], dCchanRs[-2 + nR2, 0, 0, 1, 0, 0, 2] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 2] + dPchanRs[-2 + nR2, 0, 0, 1, 0, 0, 2]]]], [[[dPchanRs[-2 + nR2, 0, 0, 2, 0, 0, 0], dPchanRs[-2 + nR2, 0, 0, 2, 0, 0, 1], dPchanRs[-2 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	cfg.dCchanRQ[R + nR2] = [[[[[[dCchanRs[-2 + nR2, 0, 0, 0, 0, 0, 0], dCchanRs[-2 + nR2, 0, 0, 0, 0, 0, 1], dCchanRs[-2 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dCchanRs[-2 + nR2, 0, 0, 1, 0, 0, 0], dCchanRs[-2 + nR2, 0, 0, 1, 0, 0, 1], dCchanRs[-2 + nR2, 0, 0, 1, 0, 0, 2] + dPchanRs[-2 + nR2, 0, 0, 1, 0, 0, 2]]]], [[[dCchanRs[-2 + nR2, 0, 0, 2, 0, 0, 0], dCchanRs[-2 + nR2, 0, 0, 2, 0, 0, 1], dCchanRs[-2 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	cfg.dDchanRQ[R + nR2] = [[[[[[dDchanRs[-2 + nR2, 0, 0, 0, 0, 0, 0], dDchanRs[-2 + nR2, 0, 0, 0, 0, 0, 1], dDchanRs[-2 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dDchanRs[-2 + nR2, 0, 0, 1, 0, 0, 0], dDchanRs[-2 + nR2, 0, 0, 1, 0, 0, 1], dDchanRs[-2 + nR2, 0, 0, 1, 0, 0, 2]]]], [[[dDchanRs[-2 + nR2, 0, 0, 2, 0, 0, 0], dDchanRs[-2 + nR2, 0, 0, 2, 0, 0, 1], dDchanRs[-2 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	R = -1
	cfg.dPchanRQ[R + nR2] = [[[[[[dCchanRs[-1 + nR2, 0, 0, 0, 0, 0, 0] + dDchanRs[nR2, 0, 0, 0, 0, 0, 0] + dPchanRs[-1 + nR2, 0, 0, 0, 0, 0, 0], dPchanRs[-1 + nR2, 0, 0, 0, 0, 0, 1], dCchanRs[-1 + nR2, 0, 0, 0, 0, 0, 2] + dDchanRs[nR2, 0, 0, 0, 0, 0, 2] + dPchanRs[-1 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dCchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0] + dPchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0], dPchanRs[-1 + nR2, 0, 0, 1, 0, 0, 1], dCchanRs[-1 + nR2, 0, 0, 1, 0, 0, 2] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 2] + dPchanRs[-1 + nR2, 0, 0, 1, 0, 0, 2]]]], [[[dPchanRs[-1 + nR2, 0, 0, 2, 0, 0, 0], dPchanRs[-1 + nR2, 0, 0, 2, 0, 0, 1], dPchanRs[-1 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	cfg.dCchanRQ[R + nR2] = [[[[[[dCchanRs[-1 + nR2, 0, 0, 0, 0, 0, 0] + dDchanRs[nR2, 0, 0, 0, 0, 0, 0] + dPchanRs[-1 + nR2, 0, 0, 0, 0, 0, 0], dCchanRs[-1 + nR2, 0, 0, 0, 0, 0, 1], dCchanRs[-1 + nR2, 0, 0, 0, 0, 0, 2] + dDchanRs[nR2, 0, 0, 0, 0, 0, 2] + dPchanRs[-1 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dCchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0] + dPchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0], dCchanRs[-1 + nR2, 0, 0, 1, 0, 0, 1] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 1], dCchanRs[-1 + nR2, 0, 0, 1, 0, 0, 2] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 2] + dPchanRs[-1 + nR2, 0, 0, 1, 0, 0, 2]]]], [[[dCchanRs[-1 + nR2, 0, 0, 2, 0, 0, 0], dCchanRs[-1 + nR2, 0, 0, 2, 0, 0, 1], dCchanRs[-1 + nR2, 0, 0, 2, 0, 0, 2] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	cfg.dDchanRQ[R + nR2] = [[[[[[dCchanRs[nR2, 0, 0, 0, 0, 0, 0] + dDchanRs[-1 + nR2, 0, 0, 0, 0, 0, 0] + dPchanRs[-1 + nR2, 0, 0, 0, 0, 0, 0], dDchanRs[-1 + nR2, 0, 0, 0, 0, 0, 1], dCchanRs[nR2, 0, 0, 0, 0, 0, 2] + dDchanRs[-1 + nR2, 0, 0, 0, 0, 0, 2] + dPchanRs[nR2, 0, 0, 0, 0, 0, 2]]]], [[[dCchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0] + dPchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0], dCchanRs[-1 + nR2, 0, 0, 1, 0, 0, 1] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 1] + dPchanRs[-2 + nR2, 0, 0, 1, 0, 0, 1], dCchanRs[-1 + nR2, 0, 0, 1, 0, 0, 2] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 2] + dPchanRs[nR2, 0, 0, 1, 0, 0, 2]]]], [[[dDchanRs[-1 + nR2, 0, 0, 2, 0, 0, 0], dDchanRs[-1 + nR2, 0, 0, 2, 0, 0, 1], dCchanRs[1 + nR2, 0, 0, 2, 0, 0, 2] + dDchanRs[-1 + nR2, 0, 0, 2, 0, 0, 2] + dPchanRs[nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	R = 0
	cfg.dPchanRQ[R + nR2] = [[[[[[dCchanRs[nR2, 0, 0, 0, 0, 0, 0] + dDchanRs[nR2, 0, 0, 0, 0, 0, 0] + dPchanRs[nR2, 0, 0, 0, 0, 0, 0], dCchanRs[nR2, 0, 0, 0, 0, 0, 1] + dDchanRs[nR2, 0, 0, 0, 0, 0, 1] + dPchanRs[nR2, 0, 0, 0, 0, 0, 1], dCchanRs[nR2, 0, 0, 0, 0, 0, 2] + dDchanRs[nR2, 0, 0, 0, 0, 0, 2] + dPchanRs[nR2, 0, 0, 0, 0, 0, 2]]]], [[[dCchanRs[nR2, 0, 0, 1, 0, 0, 0] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0] + dPchanRs[nR2, 0, 0, 1, 0, 0, 0], dCchanRs[nR2, 0, 0, 1, 0, 0, 1] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 1] + dPchanRs[nR2, 0, 0, 1, 0, 0, 1], dCchanRs[nR2, 0, 0, 1, 0, 0, 2] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 2] + dPchanRs[nR2, 0, 0, 1, 0, 0, 2]]]], [[[dCchanRs[nR2, 0, 0, 2, 0, 0, 0] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 0] + dPchanRs[nR2, 0, 0, 2, 0, 0, 0], dCchanRs[nR2, 0, 0, 2, 0, 0, 1] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 1] + dPchanRs[nR2, 0, 0, 2, 0, 0, 1], dCchanRs[nR2, 0, 0, 2, 0, 0, 2] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 2] + dPchanRs[nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	cfg.dCchanRQ[R + nR2] = [[[[[[dCchanRs[nR2, 0, 0, 0, 0, 0, 0] + dDchanRs[nR2, 0, 0, 0, 0, 0, 0] + dPchanRs[nR2, 0, 0, 0, 0, 0, 0], dCchanRs[nR2, 0, 0, 0, 0, 0, 1] + dDchanRs[nR2, 0, 0, 0, 0, 0, 1] + dPchanRs[nR2, 0, 0, 0, 0, 0, 1], dCchanRs[nR2, 0, 0, 0, 0, 0, 2] + dDchanRs[nR2, 0, 0, 0, 0, 0, 2] + dPchanRs[nR2, 0, 0, 0, 0, 0, 2]]]], [[[dCchanRs[nR2, 0, 0, 1, 0, 0, 0] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0] + dPchanRs[nR2, 0, 0, 1, 0, 0, 0], dCchanRs[nR2, 0, 0, 1, 0, 0, 1] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 1] + dPchanRs[nR2, 0, 0, 1, 0, 0, 1], dCchanRs[nR2, 0, 0, 1, 0, 0, 2] + dPchanRs[nR2, 0, 0, 1, 0, 0, 2]]]], [[[dCchanRs[nR2, 0, 0, 2, 0, 0, 0] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 0] + dPchanRs[nR2, 0, 0, 2, 0, 0, 0], dCchanRs[nR2, 0, 0, 2, 0, 0, 1] + dPchanRs[nR2, 0, 0, 2, 0, 0, 1], dCchanRs[nR2, 0, 0, 2, 0, 0, 2] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 2] + dPchanRs[nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	cfg.dDchanRQ[R + nR2] = [[[[[[dCchanRs[nR2, 0, 0, 0, 0, 0, 0] + dDchanRs[nR2, 0, 0, 0, 0, 0, 0] + dPchanRs[nR2, 0, 0, 0, 0, 0, 0], dCchanRs[nR2, 0, 0, 0, 0, 0, 1] + dDchanRs[nR2, 0, 0, 0, 0, 0, 1] + dPchanRs[-1 + nR2, 0, 0, 0, 0, 0, 1], dCchanRs[nR2, 0, 0, 0, 0, 0, 2] + dDchanRs[nR2, 0, 0, 0, 0, 0, 2] + dPchanRs[1 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dCchanRs[-1 + nR2, 0, 0, 1, 0, 0, 0] + dDchanRs[nR2, 0, 0, 1, 0, 0, 0] + dPchanRs[nR2, 0, 0, 1, 0, 0, 0], dCchanRs[-1 + nR2, 0, 0, 1, 0, 0, 1] + dDchanRs[nR2, 0, 0, 1, 0, 0, 1] + dPchanRs[-1 + nR2, 0, 0, 1, 0, 0, 1], dDchanRs[nR2, 0, 0, 1, 0, 0, 2]]]], [[[dCchanRs[1 + nR2, 0, 0, 2, 0, 0, 0] + dDchanRs[nR2, 0, 0, 2, 0, 0, 0] + dPchanRs[nR2, 0, 0, 2, 0, 0, 0], dDchanRs[nR2, 0, 0, 2, 0, 0, 1], dCchanRs[1 + nR2, 0, 0, 2, 0, 0, 2] + dDchanRs[nR2, 0, 0, 2, 0, 0, 2] + dPchanRs[1 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	R = 1
	cfg.dPchanRQ[R + nR2] = [[[[[[dCchanRs[1 + nR2, 0, 0, 0, 0, 0, 0] + dDchanRs[nR2, 0, 0, 0, 0, 0, 0] + dPchanRs[1 + nR2, 0, 0, 0, 0, 0, 0], dCchanRs[1 + nR2, 0, 0, 0, 0, 0, 1] + dDchanRs[nR2, 0, 0, 0, 0, 0, 1] + dPchanRs[1 + nR2, 0, 0, 0, 0, 0, 1], dPchanRs[1 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dPchanRs[1 + nR2, 0, 0, 1, 0, 0, 0], dPchanRs[1 + nR2, 0, 0, 1, 0, 0, 1], dPchanRs[1 + nR2, 0, 0, 1, 0, 0, 2]]]], [[[dCchanRs[1 + nR2, 0, 0, 2, 0, 0, 0] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 0] + dPchanRs[1 + nR2, 0, 0, 2, 0, 0, 0], dCchanRs[1 + nR2, 0, 0, 2, 0, 0, 1] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 1] + dPchanRs[1 + nR2, 0, 0, 2, 0, 0, 1], dPchanRs[1 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	cfg.dCchanRQ[R + nR2] = [[[[[[dCchanRs[1 + nR2, 0, 0, 0, 0, 0, 0] + dDchanRs[nR2, 0, 0, 0, 0, 0, 0] + dPchanRs[1 + nR2, 0, 0, 0, 0, 0, 0], dCchanRs[1 + nR2, 0, 0, 0, 0, 0, 1] + dDchanRs[nR2, 0, 0, 0, 0, 0, 1] + dPchanRs[1 + nR2, 0, 0, 0, 0, 0, 1], dCchanRs[1 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dCchanRs[1 + nR2, 0, 0, 1, 0, 0, 0], dCchanRs[1 + nR2, 0, 0, 1, 0, 0, 1] + dDchanRs[-1 + nR2, 0, 0, 1, 0, 0, 1], dCchanRs[1 + nR2, 0, 0, 1, 0, 0, 2]]]], [[[dCchanRs[1 + nR2, 0, 0, 2, 0, 0, 0] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 0] + dPchanRs[1 + nR2, 0, 0, 2, 0, 0, 0], dCchanRs[1 + nR2, 0, 0, 2, 0, 0, 1] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 1] + dPchanRs[1 + nR2, 0, 0, 2, 0, 0, 1], dCchanRs[1 + nR2, 0, 0, 2, 0, 0, 2] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	cfg.dDchanRQ[R + nR2] = [[[[[[dCchanRs[nR2, 0, 0, 0, 0, 0, 0] + dDchanRs[1 + nR2, 0, 0, 0, 0, 0, 0] + dPchanRs[1 + nR2, 0, 0, 0, 0, 0, 0], dCchanRs[nR2, 0, 0, 0, 0, 0, 1] + dDchanRs[1 + nR2, 0, 0, 0, 0, 0, 1] + dPchanRs[nR2, 0, 0, 0, 0, 0, 1], dDchanRs[1 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dDchanRs[1 + nR2, 0, 0, 1, 0, 0, 0], dCchanRs[-1 + nR2, 0, 0, 1, 0, 0, 1] + dDchanRs[1 + nR2, 0, 0, 1, 0, 0, 1] + dPchanRs[nR2, 0, 0, 1, 0, 0, 1], dDchanRs[1 + nR2, 0, 0, 1, 0, 0, 2]]]], [[[dCchanRs[1 + nR2, 0, 0, 2, 0, 0, 0] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 0] + dPchanRs[1 + nR2, 0, 0, 2, 0, 0, 0], dCchanRs[1 + nR2, 0, 0, 2, 0, 0, 1] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 1] + dPchanRs[nR2, 0, 0, 2, 0, 0, 1], dCchanRs[1 + nR2, 0, 0, 2, 0, 0, 2] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 2] + dPchanRs[2 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	R = 2
	cfg.dPchanRQ[R + nR2] = [[[[[[dPchanRs[2 + nR2, 0, 0, 0, 0, 0, 0], dPchanRs[2 + nR2, 0, 0, 0, 0, 0, 1], dPchanRs[2 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dPchanRs[2 + nR2, 0, 0, 1, 0, 0, 0], dPchanRs[2 + nR2, 0, 0, 1, 0, 0, 1], dPchanRs[2 + nR2, 0, 0, 1, 0, 0, 2]]]], [[[dPchanRs[2 + nR2, 0, 0, 2, 0, 0, 0], dCchanRs[2 + nR2, 0, 0, 2, 0, 0, 1] + dDchanRs[1 + nR2, 0, 0, 2, 0, 0, 1] + dPchanRs[2 + nR2, 0, 0, 2, 0, 0, 1], dPchanRs[2 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	cfg.dCchanRQ[R + nR2] = [[[[[[dCchanRs[2 + nR2, 0, 0, 0, 0, 0, 0], dCchanRs[2 + nR2, 0, 0, 0, 0, 0, 1], dCchanRs[2 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dCchanRs[2 + nR2, 0, 0, 1, 0, 0, 0], dCchanRs[2 + nR2, 0, 0, 1, 0, 0, 1], dCchanRs[2 + nR2, 0, 0, 1, 0, 0, 2]]]], [[[dCchanRs[2 + nR2, 0, 0, 2, 0, 0, 0], dCchanRs[2 + nR2, 0, 0, 2, 0, 0, 1] + dPchanRs[2 + nR2, 0, 0, 2, 0, 0, 1], dCchanRs[2 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]
	cfg.dDchanRQ[R + nR2] = [[[[[[dDchanRs[2 + nR2, 0, 0, 0, 0, 0, 0], dDchanRs[2 + nR2, 0, 0, 0, 0, 0, 1], dDchanRs[2 + nR2, 0, 0, 0, 0, 0, 2]]]], [[[dDchanRs[2 + nR2, 0, 0, 1, 0, 0, 0], dDchanRs[2 + nR2, 0, 0, 1, 0, 0, 1], dDchanRs[2 + nR2, 0, 0, 1, 0, 0, 2]]]], [[[dDchanRs[2 + nR2, 0, 0, 2, 0, 0, 0], dDchanRs[2 + nR2, 0, 0, 2, 0, 0, 1], dDchanRs[2 + nR2, 0, 0, 2, 0, 0, 2]]]]]]]


	#---------------------------------projection--------------------------


	cfg.dPchanRQ = np.reshape(cfg.dPchanRQ,(nR,NNN,NNN))
	cfg.dCchanRQ = np.reshape(cfg.dCchanRQ,(nR,NNN,NNN))
	cfg.dDchanRQ = np.reshape(cfg.dDchanRQ,(nR,NNN,NNN))

	InvFourierTrans()
	


	for qi in range(cfg.nq):
		cfg.Pchan[qi] = np.add(cfg.Pchan[qi], cfg.dPchan[qi]*dL)
		cfg.Cchan[qi] = np.add(cfg.Cchan[qi], cfg.dCchan[qi]*dL)
		cfg.Dchan[qi] = np.add(cfg.Dchan[qi], cfg.dDchan[qi]*dL)







def analysis_sm(Lc):
	pmin = 100
	cmin = 100
	dmin = 100
	qp, qc, qd = 0, 0, 0
	for q in range(cfg.nq):
		Vsc, Vsdw, Vcdw = cfg.Pchan[q], -cfg.Cchan[q], -cfg.Cchan[q] + 2 * cfg.Dchan[q]
		Pv = np.real(sorted(np.linalg.eigvals(Vsc))[0])
		Cv = np.real(sorted(np.linalg.eigvals(Vsdw))[0])
		Dv = np.real(sorted(np.linalg.eigvals(Vcdw))[0])
		if Pv < pmin:
			pmin = Pv
			qp = cfg.qrange[q]
		if Cv < cmin:
			cmin = Cv
			qc = cfg.qrange[q]
		if Dv < dmin:
			dmin = Dv
			qd = cfg.qrange[q]
	data = str([Lc, pmin, cmin, dmin])+"\n"+str(['q', qp, qc, qd])+"\n"
	data = data.replace("[","")
	data = data.replace("]","")
	data = [data, pmin, cmin, dmin]
	return data

