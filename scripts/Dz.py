# Makes a plot showing the sign of the z component of the neutrino density
# Red is parallel, blue is antiparallel, and green means at least one is zero.
import numpy as np
import matplotlib.pyplot as plt

ng = 16

r = np.genfromtxt("density_s1_g1_.txt",usecols=(0))
nr = len(r)
X, Y = np.meshgrid(r,range(ng+1))

# neutrino plot
densityZ = np.zeros((ng,nr))
for g in range(ng):
    densityZ[g,:] = np.genfromtxt("density_s1_g"+str(g+1)+"_.txt",usecols=(1)) \
                  - np.genfromtxt("density_s3_g"+str(g+1)+"_.txt",usecols=(1))
D = densityZ / abs(densityZ)
D[np.where(D != D)] = 0
plt.pcolormesh(X,Y,D)
plt.savefig("Dz_nu.pdf",bbox_inches="tight")

# antineutrino plot
plt.cla()
densityZ = np.zeros((ng,nr))
for g in range(ng):
    densityZ[g,:] = np.genfromtxt("density_s2_g"+str(g+1)+"_.txt",usecols=(1)) \
                  - np.genfromtxt("density_s3_g"+str(g+1)+"_.txt",usecols=(1))
D = densityZ / abs(densityZ)
D[np.where(D != D)] = 0
plt.pcolormesh(X,Y,D)
plt.savefig("Dz_antinu.pdf",bbox_inches="tight")
