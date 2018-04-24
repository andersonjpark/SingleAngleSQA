import numpy as np
import matplotlib.pyplot as plt

ng = 16

filename = "output/f.dat"
r = np.genfromtxt(filename,usecols=(0))
nonegroup = 2*2*2*2 #two flavors, two flavors, two states, real/imag
nonestate = 2*2*2

for m in range(2):
    plt.clf()
    for ig in range(ng):
        feeR = np.genfromtxt(filename,usecols=(1 + ig*nonegroup + m*nonestate + 0))
        fexR = np.genfromtxt(filename,usecols=(1 + ig*nonegroup + m*nonestate + 2))
        fexI = np.genfromtxt(filename,usecols=(1 + ig*nonegroup + m*nonestate + 3))
        fxxR = np.genfromtxt(filename,usecols=(1 + ig*nonegroup + m*nonestate + 6))

        plee  ,= plt.plot(r,feeR, 'b-', label=(r"$f_{ee}$" if ig==0 else ""))
        plexR ,= plt.plot(r,fexR, 'g-', label=(r"$\mathrm{Re}(f_{ex})$" if ig==0 else ""))
        plexI ,= plt.plot(r,fexI, 'r-', label=(r"$\mathrm{Im}(f_{ex})$" if ig==0 else ""))
        plxx  ,= plt.plot(r,fxxR, 'k-', label=(r"$f_{xx}$" if ig==0 else ""))

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("f_state"+str(m)+".pdf",bbox_inches="tight")

