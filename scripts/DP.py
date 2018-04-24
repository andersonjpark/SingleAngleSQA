# Makes a plot showing the relative sign of the z component of the density and potential isospin vectors.
# Red is parallel, blue is antiparallel, and green means at least one is zero.
import numpy as np
import matplotlib.pyplot as plt

ng = 16

r = np.genfromtxt("density_s1_g1_.txt",usecols=(0))
nr = len(r)
X, Y = np.meshgrid(r,range(ng+1))

# neutrino plot
densityZ = np.zeros((ng,nr))
potentialZ = np.zeros((ng,nr))
for g in range(ng):
    densityZ[g,:] = np.genfromtxt("density_s1_g"+str(g+1)+"_.txt",usecols=(1)) \
                  - np.genfromtxt("density_s3_g"+str(g+1)+"_.txt",usecols=(1))
    potentialZ[g,:] = np.genfromtxt("potential_s1_g"+str(g+1)+"_.txt",usecols=(1)) \
                    - np.genfromtxt("potential_s3_g"+str(g+1)+"_.txt",usecols=(1))
DP = densityZ*potentialZ / abs(densityZ*potentialZ)
DP[np.where(DP != DP)] = 0
plt.pcolormesh(X,Y,DP)
plt.savefig("DP_nu.pdf",bbox_inches="tight")

# antineutrino plot
plt.cla()
densityZ = np.zeros((ng,nr))
potentialZ = np.zeros((ng,nr))
for g in range(ng):
    densityZ[g,:] = np.genfromtxt("density_s2_g"+str(g+1)+"_.txt",usecols=(1)) \
                  - np.genfromtxt("density_s3_g"+str(g+1)+"_.txt",usecols=(1))
    potentialZ[g,:] = np.genfromtxt("potential_s2_g"+str(g+1)+"_.txt",usecols=(1)) \
                    - np.genfromtxt("potential_s3_g"+str(g+1)+"_.txt",usecols=(1))
DP = densityZ*potentialZ / abs(densityZ*potentialZ)
DP[np.where(DP != DP)] = 0
plt.pcolormesh(X,Y,DP)
plt.savefig("DP_antinu.pdf",bbox_inches="tight")
