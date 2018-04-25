import h5py
import numpy as np
import matplotlib.pyplot as plt

# constants
h = 6.6260755e-27 # erg s
c = 2.99792458e10 # cm/s
MeV = 1.60218e-6 # erg


# read in the interaction rate table
nulib = h5py.File("/mnt/scratch/tables/NuLib/NuLib_SFHo.h5","r")
rhogrid = nulib["rho_points"]
tempgrid = nulib["temp_points"]
yegrid = nulib["ye_points"]
absopac = np.array(nulib["absorption_opacity"])
scatopac = np.array(nulib["scattering_opacity"])
emis = np.array(nulib["emissivities"])

# get energy grid
E = np.array(nulib["neutrino_energies"])*MeV # erg
ng = len(E)
E1 = np.array(nulib["bin_top"])*MeV # erg
E0 = np.array(nulib["bin_bottom"])*MeV # erg
dE  = [E1[i]    - E0[i]    for i in range(ng)]
dE3 = [E1[i]**3 - E0[i]**3 for i in range(ng)]
ns = 3

def locateLeft(ygrid, yval):
    return next(z for z,val in enumerate(ygrid) if val > yval)-1

def interpolate(array, rho, temp, ye):
    # condition the inputs
    rho  = min(max(rho,  rhogrid[0]), rhogrid[-1])
    temp = min(max(temp,tempgrid[0]),tempgrid[-1])
    ye   = min(max(ye,    yegrid[0]),  yegrid[-1])

    # get the left indices
    irhoL  = locateLeft(rhogrid, rho)
    itempL = locateLeft(tempgrid, temp)
    iyeL   = locateLeft(yegrid, ye)
    irhoR  = irhoL+1
    itempR = itempL+1
    iyeR   = iyeL+1

    # get the left and right coordinates
    lrho0  = np.log10( rhogrid[irhoL   ])
    lrho1  = np.log10( rhogrid[irhoL +1])
    ltemp0 = np.log10(tempgrid[itempL  ])
    ltemp1 = np.log10(tempgrid[itempL+1])
    ye0    =          yegrid[iyeL      ]
    ye1    =          yegrid[iyeL    +1]

    # get the deltas
    lrho  = np.log10(rho)
    ltemp = np.log10(temp)
    dlrhoL  = abs(lrho   - lrho0)
    dlrhoR  = abs(lrho1  - lrho)
    dltempL = abs(ltemp  - ltemp0)
    dltempR = abs(ltemp1 - ltemp)
    dyeL    = abs(ye     - ye0)
    dyeR    = abs(ye1    - ye)
    vol = abs((lrho1-lrho0) * (ltemp1-ltemp0) * (ye1-ye0))

    # get the arrays at the corners
    laLLL = np.log10(array[:,:,iyeL,itempL,irhoL])
    laLLR = np.log10(array[:,:,iyeL,itempL,irhoR])
    laLRL = np.log10(array[:,:,iyeL,itempR,irhoL])
    laLRR = np.log10(array[:,:,iyeL,itempR,irhoR])
    laRLL = np.log10(array[:,:,iyeR,itempL,irhoL])
    laRLR = np.log10(array[:,:,iyeR,itempL,irhoR])
    laRRL = np.log10(array[:,:,iyeR,itempR,irhoL])
    laRRR = np.log10(array[:,:,iyeR,itempR,irhoR])

    # perform the interpolation
    la = 0
    la += laLLL * dyeR * dltempR * dlrhoR
    la += laLLR * dyeR * dltempR * dlrhoL
    la += laLRL * dyeR * dltempL * dlrhoR
    la += laLRR * dyeR * dltempL * dlrhoL
    la += laRLL * dyeL * dltempR * dlrhoR
    la += laRLR * dyeL * dltempR * dlrhoL
    la += laRRL * dyeL * dltempL * dlrhoR
    la += laRRR * dyeL * dltempL * dlrhoL

    return 10**(la/vol)
    
# get the rho,T,Ye arrays
rlist = np.genfromtxt("rho.txt", usecols=(0))
rho  = np.genfromtxt("rho.txt", usecols=(1))
temp = np.genfromtxt("temp.txt",usecols=(1))
ye   = np.genfromtxt("Ye.txt",  usecols=(1))

e = np.ndarray((len(rlist),ng,ns))
a = np.ndarray((len(rlist),ng,ns))
s = np.ndarray((len(rlist),ng,ns))
print(np.shape(e))
for ir in range(len(rlist)):
    print(ir,rlist[ir])
    s[ir,:,:] = interpolate(scatopac, rho[ir], temp[ir], ye[ir])
    a[ir,:,:] = interpolate(absopac,  rho[ir], temp[ir], ye[ir])
    e[ir,:,:] = interpolate(emis,     rho[ir], temp[ir], ye[ir])
print(np.shape(s))

# make the plots
for species in range(ns):
    plt.clf()
    for ig in range(ng):
        plt.loglog(rlist, a[:,ig,species],'k-') # 1/cm
    plt.savefig("absopac"+str(species)+".pdf",bbox_inches="tight")

    plt.clf()
    for ig in range(ng):
        plt.loglog(rlist, s[:,ig,species],'k-') # 1/cm
    plt.savefig("scatopac"+str(species)+".pdf",bbox_inches="tight")

    plt.clf()
    for ig in range(ng):
        conversion_factor = h**3 * c**2 * dE[ig] / (E[ig] * dE3[ig]/3.) / MeV
        plt.loglog(rlist, e[:,ig,species] * conversion_factor,'k-') # 1/(cm.sr)
    plt.savefig("emis"+str(species)+".pdf",bbox_inches="tight")
