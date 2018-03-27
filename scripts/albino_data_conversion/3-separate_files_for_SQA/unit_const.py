import math
import numpy as np

####################### units ############################################
MeV = 1.60217657e-6  # erg
eV = 1.60217657e-12  # erg
M = 1e6
GeV = (1.e3) * MeV  # erg
GF = (1.166379e-5) / (GeV ** 2)  # erg^{-2}
km = 1.e5  # cm
# constants
pi = 3.1415926
hcMeV = 1.23984e-10  # MeV.cm
hc = 1.23984e-10 * MeV  # erg*cm
hbarc = hc / (2 * pi)  # erg*cm
NA = 6.0221413e23
Mp = 1.67262e-24
makeAPotential = math.sqrt(2) * GF * (hbarc ** 3)  # erg^{-2}
MeVtoGK = 11.6045221
Mp = 1.67262e-24
Mn = 1.67493e-24
c2 = 8.98755e+20
c1 = math.sqrt(c2)
Me = 9.10938e-28
########################################
cV = 1.e0
cA = 1.26e0
nucleonMassDifference = ((Mn)-(Mp))*(c2)
thresholdEnergy       = nucleonMassDifference + Me*(c2)
######################test set up ############################################
# Parameters for Hydro-simulation from Albino Perego
nEnergyBins = 8
numberofneutrinotype = 3
maxE = 37.48
minE = 2.
E = np.array([2., 3.04, 4.62, 7.02, 10.67, 16.22, 24.66, 37.48])  # original bins
dE = np.array([np.exp(np.log(2*M)+(i+0.5)*((np.log(37.48*M)-np.log(2*M))/(nEnergyBins-1)))-np.exp(np.log(2*M)+(i-0.5)*((np.log(37.48*M)-np.log(2*M))/(nEnergyBins-1))) for i in range(nEnergyBins)])
dE/=M #MeV
Cnunu = math.pow(GF*hbarc**3,2)*c1*(cV*cV + cA*cA*3.e0)/(pi*math.pow(hbarc,4))