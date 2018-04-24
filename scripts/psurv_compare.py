import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

filename = "output/2p.dat"
compare_filename = "/mnt/scratch/sqa/yonglin1_orig/"+filename
r   = np.genfromtxt(        filename,usecols=(0))
rc  = np.genfromtxt(compare_filename,usecols=(0))
rmin = np.min(r)
rmax = np.max(r)
rplot = np.linspace(rmin,rmax, num=500)

pe  = interp1d(r,  np.genfromtxt(        filename,usecols=(10)))
pa  = interp1d(r,  np.genfromtxt(        filename,usecols=(11)))
pec = interp1d(rc, np.genfromtxt(compare_filename,usecols=(10)))
pac = interp1d(rc, np.genfromtxt(compare_filename,usecols=(11)))

# set up plot
fig, (ax1, ax2) = plt.subplots(2, 1)
plt.subplots_adjust(wspace=0, hspace=0)

ple,=  ax1.plot(rplot,  pe(rplot), 'b-',  label=r"$P_{\nu_e}$")
pla,=  ax1.plot(rplot,  pa(rplot), 'r-',  label=r"$P_{\bar{\nu}_e}$")
plec,= ax1.plot(rplot, pec(rplot), 'b--', label=r"$P_{\nu_e}$ (comparison)")
plac,= ax1.plot(rplot, pac(rplot), 'r--', label=r"$P_{\bar{\nu}_e}$ (comparison)")

pled,= ax2.plot(rplot, pe(rplot)-pec(rplot), 'b-')
plad,= ax2.plot(rplot, pa(rplot)-pac(rplot), 'r-')

ax1.xaxis.set_ticklabels([])
ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig("psurv_compare.pdf",bbox_inches="tight")
