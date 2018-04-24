import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

ng = 16

filename = "output/f.dat"
compare_filename = "/mnt/scratch/sqa/yonglin1_orig/"+filename
r  = np.genfromtxt(        filename,usecols=(0))
rc = np.genfromtxt(compare_filename,usecols=(0))
nonegroup = 2*2*2*2 #two flavors, two flavors, two states, real/imag
nonestate = 2*2*2

rmin = np.min(r)
rmax = np.max(r)
rplot = np.linspace(rmin,rmax, num=500)

for m in range(2):
    for ig in range(ng):
        # set up plot
        fig, (ax1, ax2) = plt.subplots(2, 1)
        plt.subplots_adjust(wspace=0, hspace=0)

        # read in data
        feeR  = interp1d(r,  np.genfromtxt(filename,        usecols=(1 + ig*nonegroup + m*nonestate + 0)) )
        fexR  = interp1d(r,  np.genfromtxt(filename,        usecols=(1 + ig*nonegroup + m*nonestate + 2)) )
        fexI  = interp1d(r,  np.genfromtxt(filename,        usecols=(1 + ig*nonegroup + m*nonestate + 3)) )
        fxxR  = interp1d(r,  np.genfromtxt(filename,        usecols=(1 + ig*nonegroup + m*nonestate + 6)) )
        feeRc = interp1d(rc, np.genfromtxt(compare_filename,usecols=(1 + ig*nonegroup + m*nonestate + 0)) )
        fexRc = interp1d(rc, np.genfromtxt(compare_filename,usecols=(1 + ig*nonegroup + m*nonestate + 2)) )
        fexIc = interp1d(rc, np.genfromtxt(compare_filename,usecols=(1 + ig*nonegroup + m*nonestate + 3)) )
        fxxRc = interp1d(rc, np.genfromtxt(compare_filename,usecols=(1 + ig*nonegroup + m*nonestate + 6)) )

        # interpolate to same locations for plotting
        plee   ,= ax1.plot(rplot,feeR(rplot),  'b-',  label=r"$f_{ee}$")
        plexR  ,= ax1.plot(rplot,fexR(rplot),  'g-',  label=r"$\mathrm{Re}(f_{ex})$")
        plexI  ,= ax1.plot(rplot,fexI(rplot),  'r-',  label=r"$\mathrm{Im}(f_{ex})$")
        plxx   ,= ax1.plot(rplot,fxxR(rplot),  'k-',  label=r"$f_{xx}$")
        pleec  ,= ax1.plot(rplot,feeRc(rplot), 'b--', label=r"$f_{ee}$ (comparison)")
        plexRc ,= ax1.plot(rplot,fexRc(rplot), 'g--', label=r"$\mathrm{Re}(f_{ex})$ (comparison)")
        plexIc ,= ax1.plot(rplot,fexIc(rplot), 'r--', label=r"$\mathrm{Im}(f_{ex})$ (comparison)")
        plxxc  ,= ax1.plot(rplot,fxxRc(rplot), 'k--', label=r"$f_{xx}$ (comparison)")

        pleed  ,= ax2.plot(rplot,feeR(rplot)-feeRc(rplot), 'b-')
        plexRd ,= ax2.plot(rplot,fexR(rplot)-fexRc(rplot), 'g-')
        plexId ,= ax2.plot(rplot,fexI(rplot)-fexIc(rplot), 'r-')
        plxxd  ,= ax2.plot(rplot,fxxR(rplot)-fxxRc(rplot), 'k-')

        # save the plot
        ax1.xaxis.set_ticklabels([])
        ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig("fcomp_state"+str(m)+"_group"+str(ig)+".pdf",bbox_inches="tight")
        plt.clf()

        
