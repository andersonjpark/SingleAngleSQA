import numpy as np
import matplotlib.pyplot as plt

filename = "output/2p.dat"
r       = np.genfromtxt(filename,usecols=(0))/1e5
pe      = np.genfromtxt(filename,usecols=(10))
pa      = np.genfromtxt(filename,usecols=(11))
peguess = np.genfromtxt(filename,usecols=(13))
paguess = np.genfromtxt(filename,usecols=(14))

plt.ylim(0,1.1)
ple,=      plt.plot(r,pe, 'b-', label=r"$P_{\nu_e}$")
pla,=      plt.plot(r,pa, 'r-', label=r"$P_{\bar{\nu}_e}$")
pleguess,= plt.plot(r, peguess, 'b--', label=r"$P_{\nu_e}$ guess")
plaguess,= plt.plot(r, paguess, 'r--', label=r"$P_{\bar{\nu}_e}$ guess")

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.xlabel("Distance (km)")
plt.ylabel("Survival Probability")

plt.savefig("survival_probabilities.pdf",bbox_inches="tight")
