import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import h5py

origin = str(45)
trajectory = str(60)
sim1 = "/mnt/scratch/sqa/fan_"+origin+"_extradata_orig/"+trajectory+"/"
sim2 = "/mnt/scratch/sqa/fan_"+origin+"_extradata_interact/"+trajectory+"/" #"/mnt/scratch/sqa/45_60_noemabs/"#
indir = "/mnt/scratch/sqa/input_data/fan_"+origin+"/"+trajectory+"/"
nulibfile = "/mnt/scratch/sqa/input_data/NuLib_rho82_temp65_ye51_ng16_ns3_version1.0_20180427.h5"
sgn = 1
    
filename1 = sim1+"output/2p.dat"
filename2 = sim2+"output/2p.dat"
pfilenamebase = "output/out.dat"
r1 = np.genfromtxt(filename1,usecols=(0))/1e5
r2 = np.genfromtxt(filename2,usecols=(0))/1e5
rmin = min(np.min(r1), np.min(r2))
rmax = 11 #max(np.max(r1), np.max(r2))
print(rmin,rmax)

# constants
h = 6.6260755e-27 # erg s
c = 2.99792458e10 # cm/s
hbarc = h*c/(2.*np.pi)
MeV = 1.60218e-6 # erg

# vacuum mixing
ng=16
grouplist = [0,int(ng/2),ng-1] #range(ng)
dm13 = 2.43e-3 # eV^2
thetaV = .15
eV = 1.60218e-12 # erg
lEtop = np.log(37.48 * 1e6) # eV
lEbottom = np.log(2. * 1e6) # eV
dlE = (lEtop-lEbottom)/(ng-1)
E =  np.array([np.exp(lEbottom + i*dlE) for i in range(ng)])
Vvac = dm13 / (2.*E) * np.cos(2.*thetaV) * eV / hbarc

# self-interaction, matter potentials, probabilities
Ve1 = np.genfromtxt(filename1,usecols=(1 ))[:len(r1)]/hbarc
VSI1 = (np.genfromtxt(filename1,usecols=(2 ))[:len(r1)] \
       - np.genfromtxt(filename1,usecols=(3 ))[:len(r1)])/hbarc # already includes anti-matter subtraction
#       - np.genfromtxt(filename1,usecols=(4 )) \
#       + np.genfromtxt(filename1,usecols=(5 ))
pe1 = np.genfromtxt(filename1,usecols=(9))[:len(r1)]/hbarc
pa1 = np.genfromtxt(filename1,usecols=(10))[:len(r1)]/hbarc

Ve2 = np.genfromtxt(filename2,usecols=(1 ))[:len(r2)]/hbarc
VSI2 = (np.genfromtxt(filename2,usecols=(2 ))[:len(r2)] \
       - np.genfromtxt(filename2,usecols=(3 ))[:len(r2)])/hbarc # already includes anti-matter subtraction
#       - np.genfromtxt(filename2,usecols=(4 )) \
#       + np.genfromtxt(filename2,usecols=(5 ))
pe2 = np.genfromtxt(filename2,usecols=(9))[:len(r2)]/hbarc
pa2 = np.genfromtxt(filename2,usecols=(10))[:len(r2)]/hbarc


# set up plot
fig, ax = plt.subplots(6, 2, squeeze=False)
fig.set_figheight(20)
fig.set_figwidth(15)
plt.subplots_adjust(wspace=0, hspace=0)

#################
# DPHI_DR PLOTS #
#################
filename1 = sim1+"output/dangledr.dat"
filename2 = sim2+"output/dangledr.dat"
rplot1 = np.genfromtxt(filename1,usecols=(0))/1e5
rplot2 = np.genfromtxt(filename2,usecols=(0))/1e5
whereinrange1 = np.where((rplot1>=rmin) & (rplot1<=rmax))
whereinrange2 = np.where((rplot2>=rmin) & (rplot2<=rmax))
phimin = 1e10
phimax = -1e10
thetamin = 1e10
thetamax = -1e10
nm = 2
def get_pos_min(array):
    wherepos = np.where(array>0)
    if len(wherepos[0])>0:
        return np.min(array[wherepos])
    else:
        return 1e100
        
for ig in grouplist:
    for im in range(nm):
        dthetadr_osc1 = sgn*np.genfromtxt(filename1, usecols=(1 + ig + im*ng + 0*nm*ng))[whereinrange1]
        dphidr_osc1   = sgn*np.genfromtxt(filename1, usecols=(1 + ig + im*ng + 1*nm*ng))[whereinrange1]
        dthetadr_int1 = sgn*np.genfromtxt(filename1, usecols=(1 + ig + im*ng + 2*nm*ng))[whereinrange1]
        dphidr_int1   = sgn*np.genfromtxt(filename1, usecols=(1 + ig + im*ng + 3*nm*ng))[whereinrange1]
        dthetadr_osc2 = sgn*np.genfromtxt(filename2, usecols=(1 + ig + im*ng + 0*nm*ng))[whereinrange2]
        dphidr_osc2   = sgn*np.genfromtxt(filename2, usecols=(1 + ig + im*ng + 1*nm*ng))[whereinrange2]
        dthetadr_int2 = sgn*np.genfromtxt(filename2, usecols=(1 + ig + im*ng + 2*nm*ng))[whereinrange2]
        dphidr_int2   = sgn*np.genfromtxt(filename2, usecols=(1 + ig + im*ng + 3*nm*ng))[whereinrange2]
        ax[3][im].semilogy(rplot1[whereinrange1], dthetadr_osc1,'--' ,color=plt.cm.winter(ig/(ng)))
        ax[3][im].semilogy(rplot1[whereinrange1], dthetadr_int1,'-',color=plt.cm.winter(ig/(ng)))
        ax[3][im].semilogy(rplot2[whereinrange2], dthetadr_osc2,'--' ,color=plt.cm.hot(ig/(ng*1.5)))
        ax[3][im].semilogy(rplot2[whereinrange2], dthetadr_int2,'-',color=plt.cm.hot(ig/(ng*1.5)))
        ax[4][im].semilogy(rplot1[whereinrange1], dphidr_osc1,  '--' ,color=plt.cm.winter(ig/(ng)))
        ax[4][im].semilogy(rplot1[whereinrange1], dphidr_int1,  '-',color=plt.cm.winter(ig/(ng)))
        ax[4][im].semilogy(rplot2[whereinrange2], dphidr_osc2,  '--' ,color=plt.cm.hot(ig/(ng*1.5)))
        ax[4][im].semilogy(rplot2[whereinrange2], dphidr_int2,  '-',color=plt.cm.hot(ig/(ng*1.5)))
        phimin = min(phimin,get_pos_min(dphidr_osc1))
        phimin = min(phimin,get_pos_min(dphidr_int1))
        phimin = min(phimin,get_pos_min(dphidr_osc2))
        phimin = min(phimin,get_pos_min(dphidr_int2))
        phimax = max(phimax,np.max(dphidr_osc1))
        phimax = max(phimax,np.max(dphidr_int1))
        phimax = max(phimax,np.max(dphidr_osc2))
        phimax = max(phimax,np.max(dphidr_int2))

        thetamin = min(thetamin,get_pos_min(dthetadr_osc1))
        thetamin = min(thetamin,get_pos_min(dthetadr_int1))
        thetamin = min(thetamin,get_pos_min(dthetadr_osc2))
        thetamin = min(thetamin,get_pos_min(dthetadr_int2))
        thetamax = max(thetamax,np.max(dthetadr_osc1))
        thetamax = max(thetamax,np.max(dthetadr_int1))
        thetamax = max(thetamax,np.max(dthetadr_osc2))
        thetamax = max(thetamax,np.max(dthetadr_int2))
        
##############################
# SURVIVAL PROBABILITY PLOTS #
##############################
def rotate_real(R,I,phi):
    return R*np.cos(phi) - I*np.sin(phi)
def rotate_imag(R,I,phi):
    return R*np.sin(phi) + I*np.cos(phi)
def norm2(R,I):
    return R*R + I*I

filename1 = sim1+pfilenamebase
filename2 = sim2+pfilenamebase
#rplot1 = np.genfromtxt(filename1,usecols=(0)) / 1e5
#rplot2 = np.genfromtxt(filename2,usecols=(0)) / 1e5
pmin = [1,1,1]
pmax = [0,0,0]
for ig in grouplist:
    for im in range(2):

        SeeR = np.genfromtxt(filename1,usecols=(1 + 8*im + 8*2*ig))
        SeeI = np.genfromtxt(filename1,usecols=(2 + 8*im + 8*2*ig))
        phi = np.arctan2(SeeI,SeeR)
        SexR = np.genfromtxt(filename1,usecols=(3 + 8*im + 8*2*ig))
        SexI = np.genfromtxt(filename1,usecols=(4 + 8*im + 8*2*ig))
        SxxR = np.genfromtxt(filename1,usecols=(7 + 8*im + 8*2*ig))
        SxxI = np.genfromtxt(filename1,usecols=(8 + 8*im + 8*2*ig))
        Pee = norm2(SeeR,SeeI)
        Pex = norm2(SexR,SexI)
        Pxx = norm2(SxxR,SxxI)

        ax[0][im].plot(rplot1, Pee, color=plt.cm.winter(ig/(ng)))
        ax[1][im].plot(rplot1, Pex, color=plt.cm.winter(ig/(ng)))
        ax[2][im].plot(rplot1, Pxx, color=plt.cm.winter(ig/(ng)))
        pmin[0] = min(pmin[0], np.min(Pee[whereinrange1]))
        pmin[1] = min(pmin[1], np.min(Pex[whereinrange1]))
        pmin[2] = min(pmin[2], np.min(Pxx[whereinrange1]))
        pmax[0] = max(pmax[0], np.max(Pee[whereinrange1]))
        pmax[1] = max(pmax[1], np.max(Pex[whereinrange1]))
        pmax[2] = max(pmax[2], np.max(Pxx[whereinrange1]))
        
        SeeR = np.genfromtxt(filename2,usecols=(1 + 8*im + 8*2*ig))
        SeeI = np.genfromtxt(filename2,usecols=(2 + 8*im + 8*2*ig))
        phi = np.arctan2(SeeI,SeeR)
        SexR = np.genfromtxt(filename2,usecols=(3 + 8*im + 8*2*ig))
        SexI = np.genfromtxt(filename2,usecols=(4 + 8*im + 8*2*ig))
        SxxR = np.genfromtxt(filename2,usecols=(7 + 8*im + 8*2*ig))
        SxxI = np.genfromtxt(filename2,usecols=(8 + 8*im + 8*2*ig))
        Pee = norm2(SeeR,SeeI)
        Pex = norm2(SexR,SexI)
        Pxx = norm2(SxxR,SxxI)

        ax[0][im].plot(rplot2, norm2(SeeR,SeeI), color=plt.cm.hot(ig/(ng*1.5)))
        ax[1][im].plot(rplot2, norm2(SexR,SexI), color=plt.cm.hot(ig/(ng*1.5)))
        ax[2][im].plot(rplot2, norm2(SxxR,SxxI), color=plt.cm.hot(ig/(ng*1.5)))
        pmin[0] = min(pmin[0], np.min(Pee[whereinrange2]))
        pmin[1] = min(pmin[1], np.min(Pex[whereinrange2]))
        pmin[2] = min(pmin[2], np.min(Pxx[whereinrange2]))
        pmax[0] = max(pmax[0], np.max(Pee[whereinrange2]))
        pmax[1] = max(pmax[1], np.max(Pex[whereinrange2]))
        pmax[2] = max(pmax[2], np.max(Pxx[whereinrange2]))
    print("Done with ig="+str(ig))

###################
# POTENTIAL PLOTS #
###################
if(any(VSI1>0)):
    pVSI1, = ax[5][0].semilogy(r1,VSI1,color="blue", label=r"$V_{\mathrm{SI},1}$")
if(any(VSI2>0)):
    pVSI2, = ax[5][0].semilogy(r2,VSI2,color="red" , label=r"$V_{\mathrm{SI},2}$")
if(any(-VSI1>0)):
    ax[5][0].semilogy(r1,-VSI1,':',color="blue", label=r"$-V_{\mathrm{SI},1}$")
if(any(-VSI2>0)):
    ax[5][0].semilogy(r2,-VSI2,':',color="red", label=r"$-V_{\mathrm{SI},2}$" )

if(r1[-1]>=r2[-1]):
    rplot = r1
    Veplot = Ve1
else:
    rplot=r2
    Veplot = Ve2

pVe,  = ax[5][0].semilogy(rplot, Veplot ,    color="black", label=r"$V_{e}$")
pVvac, = ax[5][0].semilogy(rplot, np.ones(len(rplot))*Vvac[grouplist[0]], color="gray", label=r"$V_\mathrm{vac}$")
for ig in grouplist[1:]:
    ax[5][0].semilogy(rplot, np.ones(len(rplot))*Vvac[ig], color="gray")

potential_max = np.max(Veplot)
potential_max = max(potential_max, np.max(np.abs(VSI1)))
potential_max = max(potential_max, np.max(np.abs(VSI2)))
potential_max = max(potential_max, np.max(Vvac))

potential_min = np.min(Veplot)
potential_min = min(potential_min, np.min(np.abs(VSI1)))
potential_min = min(potential_min, np.min(np.abs(VSI2)))
potential_min = min(potential_min, np.min(Vvac))

##########################
# INTERACTION RATE PLOTS #
##########################
print("Making interaction rate plots")
# read in the interaction rate table
nulib = h5py.File(nulibfile,"r")
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
rlist = np.genfromtxt(indir+"rho.txt", usecols=(0))/1e5
rho  = np.genfromtxt(indir+"rho.txt", usecols=(1))
temp = np.genfromtxt(indir+"temp.txt",usecols=(1))
ye   = np.genfromtxt(indir+"Ye.txt",  usecols=(1))

e = np.ndarray((len(rlist),ng,ns))
a = np.ndarray((len(rlist),ng,ns))
s = np.ndarray((len(rlist),ng,ns))
print(np.shape(e))
for ir in range(len(rlist)):
    s[ir,:,:] = interpolate(scatopac, rho[ir], temp[ir], ye[ir])
    a[ir,:,:] = interpolate(absopac,  rho[ir], temp[ir], ye[ir])
    e[ir,:,:] = interpolate(emis,     rho[ir], temp[ir], ye[ir])
print(np.shape(s))

# make the plots
for species in range(ns):
    if species==0:
        col = "blue"
        label = r"$\nu_e$"
    if species==1:
        col = "red"
        label = r"$\bar{\nu}_e$"
    if species==2:
        col = "black"
        label = r"$\nu_x$"
    for ig in grouplist:
        conversion_factor = h**3 * c**2 * dE[ig] / (E[ig] * dE3[ig]/3.) / MeV
        plotval = (e[:,ig,species] * conversion_factor) + a[:,ig,species] + s[:,ig,species]
        potential_max = max(potential_max,np.max(plotval))
        potential_min = min(potential_min,np.min(plotval))
        if ig==0:
            ax[5][1].semilogy(rlist, plotval,'-',color=col, label=label) # 1/cm
        else:
            ax[5][1].semilogy(rlist, plotval,'-',color=col) # 1/cm
        #ax[5][1].semilogy(rlist, s[:,ig,species],'--',color=col) # 1/cm
        #ax[5][1].semilogy(rlist, e[:,ig,species] * conversion_factor,':',color=col) # 1/(cm.sr)

################
# PLOT OPTIONS #
################
for ix in range(2):
    for iy in range(6):
        ax[iy][ix].set_xlim([rmin,rmax])
        ax[iy][ix].grid()

    for iy in range(3):
        ax[iy][ix].set_ylim([pmin[iy],pmax[iy]])
        ax[iy][ix].xaxis.set_ticklabels([])

    ax[3][ix].set_ylim([thetamin/2., thetamax*2.])
    ax[4][ix].set_ylim([phimin/2., phimax*2.])
    ax[5][ix].set_ylim([potential_min/2., potential_max*2.])
        
for iy in range(6):
    ax[iy][1].yaxis.set_label_position("right")
    ax[iy][1].yaxis.tick_right()
ax[5][0].legend(bbox_to_anchor=(-0.4, 1), loc=2, borderaxespad=0.)
ax[5][1].legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0.)
ax[0][0].set_title(sim1, color="blue", fontweight="bold")
ax[0][1].set_title(sim2, color="red", fontweight="bold")
ax[0][0].set_ylabel(r"$P_{ee}$")
ax[1][0].set_ylabel(r"$P_{ex}$")
ax[2][0].set_ylabel(r"$P_{xx}$")
ax[3][0].set_ylabel(r"$d\theta/dr$ (rad/cm)")
ax[4][0].set_ylabel(r"$d\phi/dr$ (rad/cm)")
ax[5][0].set_ylabel("Potential (cm$^{-1}$)")
ax[0][1].set_ylabel(r"$\bar{P}_{ee}$")
ax[1][1].set_ylabel(r"$\bar{P}_{ex}$")
ax[2][1].set_ylabel(r"$\bar{P}_{xx}$")
ax[3][1].set_ylabel(r"$d\theta/dr$ (rad/cm)")
ax[4][1].set_ylabel(r"$d\phi/dr$ (rad/cm)")
ax[5][1].set_ylabel("Interaction Rate (cm$^{-1}$)")
ax[5][0].set_xlabel("Distance (km)")
ax[5][1].set_xlabel("Distance (km)")
ax[0][0].xaxis.set_tick_params(labeltop='on')
ax[0][1].xaxis.set_tick_params(labeltop='on')

plt.savefig("psurv_and_potential_compare_"+origin+"_"+trajectory+".pdf",bbox_inches="tight")
