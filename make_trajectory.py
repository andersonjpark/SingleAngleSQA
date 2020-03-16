import numpy as np
import h5py

NE = 25
rho = 1e12 # g/ccm
Ye = 0.3
T = 10 # MeV

tend = 2e-5 # s

##############################################

npoints = 2 # don't change - just need the two points for an isotropic calculation
nflavors = 4
clight = 2.99792458e10 # cm/s
MeV = 1.60218e-6 # erg

Ecom = 4.*np.arange(50)+2
Etop = Ecom + 2

f = h5py.File("background.h5",'w')
f["Ecom(erg)"]    = Ecom*MeV
f["Etopcom(erg)"] = Etop*MeV
f["Elab_Elab0"] = np.ones(npoints)
f["Ecom_Elab"]  = np.ones(npoints)
f["T(MeV)"]     = np.ones(npoints) * T
f["Ye"]         = np.ones(npoints) * Ye
f["rho(g|ccm)"] = np.ones(npoints) * rho
f["ct(cm)"]     = np.array([0,tend*clight])
f["x(cm)"]      = np.array([0,tend*clight])
f["Fdens(1|ccm)"] = np.zeros((nflavors,NE,npoints))
f["Ndens(1|ccm)"] = np.zeros((nflavors,NE,npoints))
f["Pdens(1|ccm)"] = np.zeros((nflavors,NE,npoints))


f.close()
