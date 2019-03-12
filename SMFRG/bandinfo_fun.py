#Set up the band  Ek, vK, and momentum
import config as cfg 
import numpy as np 
from config import fm,fmc
import os
from config import fourierscheme

npatch = cfg.npatch	#number of patch
np4 = int(npatch/4)
np2 = int(npatch/2)
kf1, kf2 = -0.7 * cfg.pi, -0.3 * cfg.pi
tt = 2. * np.cos(kf2)
kspan = .1 * cfg.pi

NNN = cfg.NNN

nR = cfg.nR
ns = cfg.nsublattice
nf = cfg.nfactor

mu = 0

#Gf(k,i,L)[a][b]  k is the momentum, i is the sign of the frequency
# L is the frequency, a,b are the sublattice index.
def Gf(k,i,L):
	return [[1/(1j*i*L-mu+2*np.cos(k))]]

def Gf2d(k,i,L):
	return [[(-1j*i*L-2*np.cos(k))/(i**2*L**2+tt**2-4*1j*i*L*np.cos(k)-4*np.cos(k)**2),\
	tt/(i**2*L**2+tt**2-4*1j*i*L*np.cos(k)-4*np.cos(k)**2)],\
	[tt/(i**2*L**2+tt**2-4*1j*i*L*np.cos(k)-4*np.cos(k)**2),\
	(-1j*i*L-2*np.cos(k))/(i**2*L**2+tt**2-4*1j*i*L*np.cos(k)-4*np.cos(k)**2)]]



def check():
	nR = cfg.nR
	ns = cfg.nsublattice
	nf = cfg.nfactor
	cfg.PchanR = np.reshape(cfg.PchanR,(nR,ns,ns,nf,ns,ns,nf))
	cfg.CchanR = np.reshape(cfg.CchanR,(nR,ns,ns,nf,ns,ns,nf))
	cfg.DchanR = np.reshape(cfg.DchanR,(nR,ns,ns,nf,ns,ns,nf))
	
	PchanR, CchanR, DchanR = cfg.PchanR, cfg.CchanR, cfg.DchanR
	
	for ii in range(cfg.nR):
		cfg.PchanR[ii] = [[[[[[PchanR[ii, 0, 0, 0, 0, 0, 0], PchanR[ii, 0, 0, 0, 0, 0, 1], PchanR[ii, 0, 0, 0, 0, 0, 2], 0], [0, 0, 0, PchanR[ii, 0, 0, 0, 0, 1, 3]]], [[0, 0, 0, PchanR[ii, 0, 0, 0, 1, 0, 3]], [PchanR[ii, 0, 0, 0, 1, 1, 0], PchanR[ii, 0, 0, 0, 1, 1, 1], PchanR[ii, 0, 0, 0, 1, 1, 2], 0]]], [[[PchanR[ii, 0, 0, 1, 0, 0, 0], PchanR[ii, 0, 0, 1, 0, 0, 1], PchanR[ii, 0, 0, 1, 0, 0, 2], 0], [0, 0, 0, PchanR[ii, 0, 0, 1, 0, 1, 3]]], [[0, 0, 0, PchanR[ii, 0, 0, 1, 1, 0, 3]], [PchanR[ii, 0, 0, 1, 1, 1, 0], PchanR[ii, 0, 0, 1, 1, 1, 1], PchanR[ii, 0, 0, 1, 1, 1, 2], 0]]], [[[PchanR[ii, 0, 0, 2, 0, 0, 0], PchanR[ii, 0, 0, 2, 0, 0, 1], PchanR[ii, 0, 0, 2, 0, 0, 2], 0], [0, 0, 0, PchanR[ii, 0, 0, 2, 0, 1, 3]]], [[0, 0, 0, PchanR[ii, 0, 0, 2, 1, 0, 3]], [PchanR[ii, 0, 0, 2, 1, 1, 0], PchanR[ii, 0, 0, 2, 1, 1, 1], PchanR[ii, 0, 0, 2, 1, 1, 2], 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]]], [[[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[PchanR[ii, 0, 1, 3, 0, 0, 0], PchanR[ii, 0, 1, 3, 0, 0, 1], PchanR[ii, 0, 1, 3, 0, 0, 2], 0], [0, 0, 0, PchanR[ii, 0, 1, 3, 0, 1, 3]]], [[0, 0, 0, PchanR[ii, 0, 1, 3, 1, 0, 3]], [PchanR[ii, 0, 1, 3, 1, 1, 0], PchanR[ii, 0, 1, 3, 1, 1, 1], PchanR[ii, 0, 1, 3, 1, 1, 2], 0]]]]], [[[[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[PchanR[ii, 1, 0, 3, 0, 0, 0], PchanR[ii, 1, 0, 3, 0, 0, 1], PchanR[ii, 1, 0, 3, 0, 0, 2], 0], [0, 0, 0, PchanR[ii, 1, 0, 3, 0, 1, 3]]], [[0, 0, 0, PchanR[ii, 1, 0, 3, 1, 0, 3]], [PchanR[ii, 1, 0, 3, 1, 1, 0], PchanR[ii, 1, 0, 3, 1, 1, 1], PchanR[ii, 1, 0, 3, 1, 1, 2], 0]]]], [[[[PchanR[ii, 1, 1, 0, 0, 0, 0], PchanR[ii, 1, 1, 0, 0, 0, 1], PchanR[ii, 1, 1, 0, 0, 0, 2], 0], [0, 0, 0, PchanR[ii, 1, 1, 0, 0, 1, 3]]], [[0, 0, 0, PchanR[ii, 1, 1, 0, 1, 0, 3]], [PchanR[ii, 1, 1, 0, 1, 1, 0], PchanR[ii, 1, 1, 0, 1, 1, 1], PchanR[ii, 1, 1, 0, 1, 1, 2], 0]]], [[[PchanR[ii, 1, 1, 1, 0, 0, 0], PchanR[ii, 1, 1, 1, 0, 0, 1], PchanR[ii, 1, 1, 1, 0, 0, 2], 0], [0, 0, 0, PchanR[ii, 1, 1, 1, 0, 1, 3]]], [[0, 0, 0, PchanR[ii, 1, 1, 1, 1, 0, 3]], [PchanR[ii, 1, 1, 1, 1, 1, 0], PchanR[ii, 1, 1, 1, 1, 1, 1], PchanR[ii, 1, 1, 1, 1, 1, 2], 0]]], [[[PchanR[ii, 1, 1, 2, 0, 0, 0], PchanR[ii, 1, 1, 2, 0, 0, 1], PchanR[ii, 1, 1, 2, 0, 0, 2], 0], [0, 0, 0, PchanR[ii, 1, 1, 2, 0, 1, 3]]], [[0, 0, 0, PchanR[ii, 1, 1, 2, 1, 0, 3]], [PchanR[ii, 1, 1, 2, 1, 1, 0], PchanR[ii, 1, 1, 2, 1, 1, 1], PchanR[ii, 1, 1, 2, 1, 1, 2], 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]]]]]
		cfg.CchanR[ii] = [[[[[[CchanR[ii, 0, 0, 0, 0, 0, 0], CchanR[ii, 0, 0, 0, 0, 0, 1], CchanR[ii, 0, 0, 0, 0, 0, 2], 0], [0, 0, 0, CchanR[ii, 0, 0, 0, 0, 1, 3]]], [[0, 0, 0, CchanR[ii, 0, 0, 0, 1, 0, 3]], [CchanR[ii, 0, 0, 0, 1, 1, 0], CchanR[ii, 0, 0, 0, 1, 1, 1], CchanR[ii, 0, 0, 0, 1, 1, 2], 0]]], [[[CchanR[ii, 0, 0, 1, 0, 0, 0], CchanR[ii, 0, 0, 1, 0, 0, 1], CchanR[ii, 0, 0, 1, 0, 0, 2], 0], [0, 0, 0, CchanR[ii, 0, 0, 1, 0, 1, 3]]], [[0, 0, 0, CchanR[ii, 0, 0, 1, 1, 0, 3]], [CchanR[ii, 0, 0, 1, 1, 1, 0], CchanR[ii, 0, 0, 1, 1, 1, 1], CchanR[ii, 0, 0, 1, 1, 1, 2], 0]]], [[[CchanR[ii, 0, 0, 2, 0, 0, 0], CchanR[ii, 0, 0, 2, 0, 0, 1], CchanR[ii, 0, 0, 2, 0, 0, 2], 0], [0, 0, 0, CchanR[ii, 0, 0, 2, 0, 1, 3]]], [[0, 0, 0, CchanR[ii, 0, 0, 2, 1, 0, 3]], [CchanR[ii, 0, 0, 2, 1, 1, 0], CchanR[ii, 0, 0, 2, 1, 1, 1], CchanR[ii, 0, 0, 2, 1, 1, 2], 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]]], [[[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[CchanR[ii, 0, 1, 3, 0, 0, 0], CchanR[ii, 0, 1, 3, 0, 0, 1], CchanR[ii, 0, 1, 3, 0, 0, 2], 0], [0, 0, 0, CchanR[ii, 0, 1, 3, 0, 1, 3]]], [[0, 0, 0, CchanR[ii, 0, 1, 3, 1, 0, 3]], [CchanR[ii, 0, 1, 3, 1, 1, 0], CchanR[ii, 0, 1, 3, 1, 1, 1], CchanR[ii, 0, 1, 3, 1, 1, 2], 0]]]]], [[[[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[CchanR[ii, 1, 0, 3, 0, 0, 0], CchanR[ii, 1, 0, 3, 0, 0, 1], CchanR[ii, 1, 0, 3, 0, 0, 2], 0], [0, 0, 0, CchanR[ii, 1, 0, 3, 0, 1, 3]]], [[0, 0, 0, CchanR[ii, 1, 0, 3, 1, 0, 3]], [CchanR[ii, 1, 0, 3, 1, 1, 0], CchanR[ii, 1, 0, 3, 1, 1, 1], CchanR[ii, 1, 0, 3, 1, 1, 2], 0]]]], [[[[CchanR[ii, 1, 1, 0, 0, 0, 0], CchanR[ii, 1, 1, 0, 0, 0, 1], CchanR[ii, 1, 1, 0, 0, 0, 2], 0], [0, 0, 0, CchanR[ii, 1, 1, 0, 0, 1, 3]]], [[0, 0, 0, CchanR[ii, 1, 1, 0, 1, 0, 3]], [CchanR[ii, 1, 1, 0, 1, 1, 0], CchanR[ii, 1, 1, 0, 1, 1, 1], CchanR[ii, 1, 1, 0, 1, 1, 2], 0]]], [[[CchanR[ii, 1, 1, 1, 0, 0, 0], CchanR[ii, 1, 1, 1, 0, 0, 1], CchanR[ii, 1, 1, 1, 0, 0, 2], 0], [0, 0, 0, CchanR[ii, 1, 1, 1, 0, 1, 3]]], [[0, 0, 0, CchanR[ii, 1, 1, 1, 1, 0, 3]], [CchanR[ii, 1, 1, 1, 1, 1, 0], CchanR[ii, 1, 1, 1, 1, 1, 1], CchanR[ii, 1, 1, 1, 1, 1, 2], 0]]], [[[CchanR[ii, 1, 1, 2, 0, 0, 0], CchanR[ii, 1, 1, 2, 0, 0, 1], CchanR[ii, 1, 1, 2, 0, 0, 2], 0], [0, 0, 0, CchanR[ii, 1, 1, 2, 0, 1, 3]]], [[0, 0, 0, CchanR[ii, 1, 1, 2, 1, 0, 3]], [CchanR[ii, 1, 1, 2, 1, 1, 0], CchanR[ii, 1, 1, 2, 1, 1, 1], CchanR[ii, 1, 1, 2, 1, 1, 2], 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]]]]]
		cfg.DchanR[ii] = [[[[[[DchanR[ii, 0, 0, 0, 0, 0, 0], DchanR[ii, 0, 0, 0, 0, 0, 1], DchanR[ii, 0, 0, 0, 0, 0, 2], 0], [0, 0, 0, DchanR[ii, 0, 0, 0, 0, 1, 3]]], [[0, 0, 0, DchanR[ii, 0, 0, 0, 1, 0, 3]], [DchanR[ii, 0, 0, 0, 1, 1, 0], DchanR[ii, 0, 0, 0, 1, 1, 1], DchanR[ii, 0, 0, 0, 1, 1, 2], 0]]], [[[DchanR[ii, 0, 0, 1, 0, 0, 0], DchanR[ii, 0, 0, 1, 0, 0, 1], DchanR[ii, 0, 0, 1, 0, 0, 2], 0], [0, 0, 0, DchanR[ii, 0, 0, 1, 0, 1, 3]]], [[0, 0, 0, DchanR[ii, 0, 0, 1, 1, 0, 3]], [DchanR[ii, 0, 0, 1, 1, 1, 0], DchanR[ii, 0, 0, 1, 1, 1, 1], DchanR[ii, 0, 0, 1, 1, 1, 2], 0]]], [[[DchanR[ii, 0, 0, 2, 0, 0, 0], DchanR[ii, 0, 0, 2, 0, 0, 1], DchanR[ii, 0, 0, 2, 0, 0, 2], 0], [0, 0, 0, DchanR[ii, 0, 0, 2, 0, 1, 3]]], [[0, 0, 0, DchanR[ii, 0, 0, 2, 1, 0, 3]], [DchanR[ii, 0, 0, 2, 1, 1, 0], DchanR[ii, 0, 0, 2, 1, 1, 1], DchanR[ii, 0, 0, 2, 1, 1, 2], 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]]], [[[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[DchanR[ii, 0, 1, 3, 0, 0, 0], DchanR[ii, 0, 1, 3, 0, 0, 1], DchanR[ii, 0, 1, 3, 0, 0, 2], 0], [0, 0, 0, DchanR[ii, 0, 1, 3, 0, 1, 3]]], [[0, 0, 0, DchanR[ii, 0, 1, 3, 1, 0, 3]], [DchanR[ii, 0, 1, 3, 1, 1, 0], DchanR[ii, 0, 1, 3, 1, 1, 1], DchanR[ii, 0, 1, 3, 1, 1, 2], 0]]]]], [[[[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]], [[[DchanR[ii, 1, 0, 3, 0, 0, 0], DchanR[ii, 1, 0, 3, 0, 0, 1], DchanR[ii, 1, 0, 3, 0, 0, 2], 0], [0, 0, 0, DchanR[ii, 1, 0, 3, 0, 1, 3]]], [[0, 0, 0, DchanR[ii, 1, 0, 3, 1, 0, 3]], [DchanR[ii, 1, 0, 3, 1, 1, 0], DchanR[ii, 1, 0, 3, 1, 1, 1], DchanR[ii, 1, 0, 3, 1, 1, 2], 0]]]], [[[[DchanR[ii, 1, 1, 0, 0, 0, 0], DchanR[ii, 1, 1, 0, 0, 0, 1], DchanR[ii, 1, 1, 0, 0, 0, 2], 0], [0, 0, 0, DchanR[ii, 1, 1, 0, 0, 1, 3]]], [[0, 0, 0, DchanR[ii, 1, 1, 0, 1, 0, 3]], [DchanR[ii, 1, 1, 0, 1, 1, 0], DchanR[ii, 1, 1, 0, 1, 1, 1], DchanR[ii, 1, 1, 0, 1, 1, 2], 0]]], [[[DchanR[ii, 1, 1, 1, 0, 0, 0], DchanR[ii, 1, 1, 1, 0, 0, 1], DchanR[ii, 1, 1, 1, 0, 0, 2], 0], [0, 0, 0, DchanR[ii, 1, 1, 1, 0, 1, 3]]], [[0, 0, 0, DchanR[ii, 1, 1, 1, 1, 0, 3]], [DchanR[ii, 1, 1, 1, 1, 1, 0], DchanR[ii, 1, 1, 1, 1, 1, 1], DchanR[ii, 1, 1, 1, 1, 1, 2], 0]]], [[[DchanR[ii, 1, 1, 2, 0, 0, 0], DchanR[ii, 1, 1, 2, 0, 0, 1], DchanR[ii, 1, 1, 2, 0, 0, 2], 0], [0, 0, 0, DchanR[ii, 1, 1, 2, 0, 1, 3]]], [[0, 0, 0, DchanR[ii, 1, 1, 2, 1, 0, 3]], [DchanR[ii, 1, 1, 2, 1, 1, 0], DchanR[ii, 1, 1, 2, 1, 1, 1], DchanR[ii, 1, 1, 2, 1, 1, 2], 0]]], [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]]]]]

	cfg.PchanR = np.reshape(cfg.PchanR,(nR,NNN,NNN))
	cfg.CchanR = np.reshape(cfg.CchanR,(nR,NNN,NNN))
	cfg.DchanR = np.reshape(cfg.DchanR,(nR,NNN,NNN))


#define the sublattice band structure
def xi(k):
	return [-2 * np.cos(k) + tt, -2 * np.cos(k) - tt]


def InvFourierTrans_ini():
	cfg.Pchan, cfg.Cchan, cfg.Dchan = np.zeros((cfg.nq,NNN,NNN)),\
	np.zeros((cfg.nq,NNN,NNN)), np.zeros((cfg.nq,NNN,NNN))

	for m in range(cfg.NNN): 
		for n in range(cfg.NNN):
			for q in range(cfg.nq): # tabulate m,n,q
			#--------------------------------------
				for Ri in range(cfg.nR-1): #sum over R
					cfg.Pchan[q][m][n] = np.real(cfg.Pchan[q][m][n] +\
					cfg.PchanR[Ri][m][n] * np.exp(-1j*cfg.qrange[q]*cfg.R[Ri]))

					cfg.Cchan[q][m][n] = np.real(cfg.Cchan[q][m][n] +\
					cfg.CchanR[Ri][m][n] * np.exp(-1j*cfg.qrange[q]*cfg.R[Ri]))
	
					cfg.Dchan[q][m][n] = np.real(cfg.Dchan[q][m][n] +\
					cfg.DchanR[Ri][m][n] * np.exp(-1j*cfg.qrange[q]*cfg.R[Ri]))






def setbandinfo():

	#----------------------------BAND INFO-------------------------


	for i in range(np4):
		cfg.FermiM[i] = kf1 - kspan + (2*i+1) * kspan / (npatch/4)
		cfg.FermiE[i] = -2. * np.cos(cfg.FermiM[i]) - tt
		cfg.FermiV[i] = 2. * np.sin(cfg.FermiM[i])
	for i in range(np4, np2):
		i1 = np2 - i - 1
		cfg.FermiM[i] = -cfg.FermiM[i1]
		cfg.FermiE[i] = cfg.FermiE[i1]
		cfg.FermiV[i] = -cfg.FermiV[i1]
		
	dk = kf2 - kf1
	cfg.FermiM[np2:np2+np4] = cfg.FermiM[0:np4] + dk
	cfg.FermiM[np2+np4:npatch] = cfg.FermiM[np4:np2] - dk
	cfg.FermiV[np2:npatch] = cfg.FermiV[0:np2]
	cfg.FermiE[np2:npatch] = cfg.FermiE[0:np2]
	#----------------------------BAND INFO-------------------------



	#---------------------------Momentum conservation leg4---------
	for i3 in range(npatch):
		for i2 in range(npatch):
			for i1 in range(npatch):
				kv = cfg.FermiM[i1] + cfg.FermiM[i2] - cfg.FermiM[i3]
				kv = kv % cfg.Qvec
				if kv > cfg.pi:
					kv = kv - cfg.twopi
				elif kv < -cfg.pi :
					kv = kv + cfg.twopi

				for i4 in range(npatch):
					if abs((kv - cfg.FermiM[i4]) % cfg.twopi) > 1.0e-4:
						continue
					cfg.leg4[i1,i2,i3] = i4
					
	#---------------------------Momentum conservation leg4---------


	#---------------------------Assign vertex-----------------------
	U = 0.25/2;  Vnn = 0.10/2;  Jnn = 4 * ( U + Vnn )           #SC + SDW degenerate manifold, or SZH SO(5) manifold 
	#U = 0.3;  Vnn = 0;  Jnn = 0.4                              #SC + CDW  degenerate manifold
	#U = 0.3/2;  Vnn = -0.1/2;  Jnn = 0.8/2                     #SDW + CDW degenerate manifold
	Vnn = Vnn - Jnn/4
	U = 1/3

	def sgn(i):
		if i < np2:
			return 1
		else:
			return -1


	for i3 in range(npatch):
		for i2 in range(npatch):
			for i1 in range(npatch):
				i4 = cfg.leg4[i1,i2,i3]
				if i4 < 0 :
					continue
				cfg.Ga[i1,i2,i3] = cfg.Ga[i1, i2, i3] - U * ( 1 + sgn(i1) \
					* sgn(i2) * sgn(i3) * sgn(i4) ) / 4.
				cfg.Ga[i1, i2, i3] = cfg.Ga[i1, i2, i3] - Vnn * ( sgn(i2) \
				 * sgn(i3) + sgn(i1) * sgn(i4) ) / 4.       
				  #coulomb interaction on rung
				cfg.Ga[i1, i2, i3] = cfg.Ga[i1, i2, i3] + 0.5 * Jnn \
				* ( sgn(i1) * sgn(i3) + sgn(i2) * sgn(i4) ) / 4.  
				#antiferro spin exchange on rung (written as hop * hop which introduces a density-density interaction)

	#---------------------------Assign vertex-----------------------

	#-------------------------Assign channel-----------------
	#### the on-site interaction for sdw channel is -U while others are 0.
	nR = cfg.nR
	ns = cfg.nsublattice
	nf = cfg.nfactor
	NNN = cfg.NNN
	for m in range(NNN):
		for n in range(NNN):
			#for ri in range(cfg.nR):
				#cfg.PchanR[ri,m,n] = 0
				#cfg.CchanR[ri,m,n] = -U
				#cfg.DchanR[ri,m,n] = -U/2
			cfg.PchanR[int((nR-1)/2),m,n] = U
			cfg.CchanR[int((nR-1)/2),m,n] = U
			cfg.DchanR[int((nR-1)/2),m,n] = U

	
	#----------------------------------------------------
	cfg.PchanR = np.reshape(cfg.PchanR,(nR,ns,ns,nf,ns,ns,nf))
	cfg.CchanR = np.reshape(cfg.CchanR,(nR,ns,ns,nf,ns,ns,nf))
	cfg.DchanR = np.reshape(cfg.DchanR,(nR,ns,ns,nf,ns,ns,nf))
	
	PchanR, CchanR, DchanR = cfg.PchanR, cfg.CchanR, cfg.DchanR
	
	for ii in range(cfg.nR):
		cfg.PchanR[ii] = [[[[[[PchanR[ii, 0, 0, 0, 0, 0, 0], PchanR[ii, 0, 0, 0, 0, 0, 1], PchanR[ii, 0, 0, 0, 0, 0, 2]]]], [[[PchanR[ii, 0, 0, 1, 0, 0, 0], PchanR[ii, 0, 0, 1, 0, 0, 1], PchanR[ii, 0, 0, 1, 0, 0, 2]]]], [[[PchanR[ii, 0, 0, 2, 0, 0, 0], PchanR[ii, 0, 0, 2, 0, 0, 1], PchanR[ii, 0, 0, 2, 0, 0, 2]]]]]]]
		cfg.CchanR[ii] = [[[[[[CchanR[ii, 0, 0, 0, 0, 0, 0], CchanR[ii, 0, 0, 0, 0, 0, 1], CchanR[ii, 0, 0, 0, 0, 0, 2]]]], [[[CchanR[ii, 0, 0, 1, 0, 0, 0], CchanR[ii, 0, 0, 1, 0, 0, 1], CchanR[ii, 0, 0, 1, 0, 0, 2]]]], [[[CchanR[ii, 0, 0, 2, 0, 0, 0], CchanR[ii, 0, 0, 2, 0, 0, 1], CchanR[ii, 0, 0, 2, 0, 0, 2]]]]]]]
		cfg.DchanR[ii] = [[[[[[DchanR[ii, 0, 0, 0, 0, 0, 0], DchanR[ii, 0, 0, 0, 0, 0, 1], DchanR[ii, 0, 0, 0, 0, 0, 2]]]], [[[DchanR[ii, 0, 0, 1, 0, 0, 0], DchanR[ii, 0, 0, 1, 0, 0, 1], DchanR[ii, 0, 0, 1, 0, 0, 2]]]], [[[DchanR[ii, 0, 0, 2, 0, 0, 0], DchanR[ii, 0, 0, 2, 0, 0, 1], DchanR[ii, 0, 0, 2, 0, 0, 2]]]]]]]


	cfg.PchanR = np.reshape(cfg.PchanR,(nR,NNN,NNN))
	cfg.CchanR = np.reshape(cfg.CchanR,(nR,NNN,NNN))
	cfg.DchanR = np.reshape(cfg.DchanR,(nR,NNN,NNN))
	#-------------------------------------------------------

	InvFourierTrans_ini()


#	cfg.Pchan = -0.1*np.ones((cfg.nq,cfg.NNN,cfg.NNN))
#	cfg.Cchan = -0.1*np.ones((cfg.nq,cfg.NNN,cfg.NNN))
#	cfg.Dchan = -0.1*np.ones((cfg.nq,cfg.NNN,cfg.NNN))
	#-----------------------






	return "band info is OK"

