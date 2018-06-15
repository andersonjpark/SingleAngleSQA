# designed to be run in input data directory

import matplotlib.pyplot as plt
import numpy as np
import os

ng = 16
xmin = 0
xmax = 25

# vacuum mixing
dm13 = 2.43e-3 # eV^2
thetaV = .15
eV = 1.60218e-12 # erg
lEtop = np.log(37.48 * 1e6) # eV
lEbottom = np.log(2. * 1e6) # eV
dlE = (lEtop-lEbottom)/(ng-1)
E =  np.array([np.exp(lEbottom + i*dlE) for i in range(ng)])
Vvac = dm13 / (2.*E) * np.cos(2.*thetaV) * eV

# matter potential
GeV = 0.00160218 # erg
GF = 1.1663787e-5 # GeV^-2
hbar=1.05457266e-27 # erg s
c = 2.99792458e10 # cm/s
GFhbarc3 = GF * GeV**(-2) * (hbar*c)**3
mp = 1.6726219e-24 # g
rho = np.genfromtxt("rho.txt",usecols=(1))
Ye = np.genfromtxt("Ye.txt",usecols=(1))
Vmatter = np.sqrt(2) * GFhbarc3 * rho*Ye / mp

x = np.genfromtxt("potential_s1_g1_.txt",usecols=(0))/1e5
ymin = 1
ymax = 0
sumVnu = np.zeros(len(x))

for g in range(ng):
    potential = np.genfromtxt("potential_s1_g"+str(g+1)+"_.txt",usecols=(1)) \
                - np.genfromtxt("potential_s3_g"+str(g+1)+"_.txt",usecols=(1)) \
                - np.genfromtxt("potential_s2_g"+str(g+1)+"_.txt",usecols=(1)) \
                + np.genfromtxt("potential_s3_g"+str(g+1)+"_.txt",usecols=(1))
    p, =plt.semilogy(x,potential,color=plt.cm.hot(g/(ng*1.5)))
    plt.semilogy(x,-potential,":",color=plt.cm.hot(g/(ng*1.5)))
    plt.semilogy(x,np.repeat(Vvac[g],len(x)),color=plt.cm.hot(g/(ng*1.5)))
    ymin = min(ymin, np.min(abs(potential[np.where((x>=xmin) & (x<=xmax))])))# & (abs(potential)>0))])))
    ymax = max(ymax, np.max(abs(potential[np.where((x>=xmin) & (x<=xmax))])))
    sumVnu = sumVnu + potential

plt.semilogy(x,Vmatter,color="green", label="Vmatter")
plt.semilogy(x,sumVnu,color="blue", label="Vnu")
plt.semilogy(x,-sumVnu,":",color="blue", label="Vnu")
ymin = min(ymin, Vvac[ng-1])
ymax = max(ymax, np.max(Vmatter))
plt.legend()
plt.xlabel("Distance (km)")
plt.ylabel("Potential")
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.grid()
plt.title(os.getcwd())
plt.savefig("potentials.pdf")
        
