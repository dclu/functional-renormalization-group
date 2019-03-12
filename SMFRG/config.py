###The module file for all global variables  config.py
import numpy as np 

extra = 0
npatch = int((1 + 2 * extra)*4)	#number of patch

leg4 = -np.ones((npatch,npatch,npatch))		#momentum conservation
leg4 = leg4.astype(np.int)

xpp, xph = np.zeros((npatch,npatch)), np.zeros((npatch,npatch)) 	#bubble
Ga, dGa = np.zeros((npatch,npatch,npatch)), np.zeros((npatch,npatch,npatch))	#vertex function for frg
FermiE, FermiV, FermiM = np.zeros(npatch), np.zeros(npatch), np.zeros(npatch) 	#Fermi energy velocity and momentum
twopi, pi = np.arcsin(1) * 4, np.arcsin(1) * 2

Qvec = twopi



### Define for SMFRG

nq = 25
nR = 11
qrange = (np.array(range(nq))-(nq-1)/2*np.ones(nq))/((nq-1)/2)*pi


#fm(k)[a][b][m]  form factor, a,b are sublattice index,\
# m is onsite, neighbour. 
def fm(k):
	res = np.array([[[1,np.exp(-1j*k),np.exp(1j*k),0],[0,0,0,1]],\
		[[0,0,0,1],[1,np.exp(-1j*k),np.exp(1j*k),0]]])/twopi
	return res


def fmc(k):
	res = np.array([[[1,np.exp(1j*k),np.exp(-1j*k),0],[0,0,0,1]],\
		[[0,0,0,1],[1,np.exp(1j*k),np.exp(-1j*k),0]]])
	return res



nfactor = 3
nsublattice = 1


NNN = nfactor*nsublattice**2

#xpps[q,(abm),(cdn)], xphs[q,(abm),(cdn)]
xpps, xphs = np.zeros((nq,NNN,NNN)), np.zeros((nq,NNN,NNN))

#Kchan[q,(abm),(cdn)]
Pchan, Cchan, Dchan = np.zeros((nq,NNN,NNN)),\
np.zeros((nq,NNN,NNN)), np.zeros((nq,NNN,NNN))

#dKchan[q,(abm),(cdn)]
dPchan, dCchan, dDchan = np.zeros((nq,NNN,NNN)),\
np.zeros((nq,NNN,NNN)), np.zeros((nq,NNN,NNN))






###Find allowed projection,
P2C, C2P = [],[]
P2D, D2P = [],[]
C2D, D2C = [],[]
cut = 0.01

def allowed_proj():
	for R_p in range(-10,10+1):
		for rp_p in range(-1,1+1):
			for r_p in range(-1,1+1):
				for R_c in range(-10,10+1):
					for rp_c in range(-1,1+1):
						for r_c in range(-1,1+1):
							if abs(R_c + rp_c - r_p) <= cut\
							and abs(R_c - R_p) <= cut\
							and abs(R_p + rp_p - r_c) <= cut:
								P2C.append([R_p,r_p,rp_p,R_c])
								C2P.append([R_c,r_c,rp_c,R_p])



	for R_p in range(-10,10+1):
		for rp_p in range(-1,1+1):
			for r_p in range(-1,1+1):
				for R_d in range(-10,10+1):
					for rp_d in range(-1,1+1):
						for r_d in range(-1,1+1):
							if abs(R_d + rp_d - r_p) <= cut\
							and abs(R_p - r_d) <= cut\
							and abs(R_p + rp_p - R_d) <= cut:
								P2D.append([R_p,r_p,rp_p,R_d])
								D2P.append([R_d,r_d,rp_d,R_p])


	for R_c in range(-10,10+1):
		for rp_c in range(-1,1+1):
			for r_c in range(-1,1+1):
				for R_d in range(-10,10+1):
					for rp_d in range(-1,1+1):
						for r_d in range(-1,1+1):
							if abs(R_d + rp_d - R_c - rp_c) <= cut\
							and abs(R_d - r_c) <= cut\
							and abs(r_d - R_c) <= cut:
								C2D.append([R_c,r_c,rp_c,R_d])
								D2C.append([R_d,r_d,rp_d,R_c])




#   0 -->------>--- R
#           |
#           |
#   r -->------>--- R + rp


nr = 3
R = np.array(range(nR))-(nR-1)/2*np.ones(nR)
r = np.array([-1, 0, 1])
rp = np.array([-1, 0, 1])
R = R.astype(int)
r = r.astype(int)
rp = rp.astype(int)

### partial increasement
## dPchanR[nR,(abm),(cdn)]

dPchanR, dCchanR,dDchanR = np.zeros((nR,NNN,NNN),dtype=complex),\
np.zeros((nR,NNN,NNN),dtype=complex), np.zeros((nR,NNN,NNN),dtype=complex)



### total increasement
dPchanRQ, dCchanRQ,dDchanRQ = np.zeros((nR,NNN,NNN),dtype=complex),\
np.zeros((nR,NNN,NNN),dtype=complex), np.zeros((nR,NNN,NNN),dtype=complex)



PchanR, CchanR,DchanR = np.zeros((nR,NNN,NNN),dtype=complex),\
np.zeros((nR,NNN,NNN),dtype=complex), np.zeros((nR,NNN,NNN),dtype=complex)




bubblescheme = 'integration'
fourierscheme = 'ldc'