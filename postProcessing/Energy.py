import math
import numpy as np
import matplotlib.pylab as plt

U=np.loadtxt("U_1.dat")
#K=np.loadtxt("K_1.dat")
labels=["Uion","Ugas","Uvap","Ugi","Ugg","Uvg","Uvi","Uvv"]

for i in np.arange(np.size(U[0])-1):
	plt.plot(U.T[0]*1e-6,U.T[i+1],label=labels[i])

plt.plot(U.T[0]*1e-6,np.sum(U.T[1:],axis=0),label="Total")
plt.xlabel("Time [ns]")
plt.ylabel("Energy [kcal/mol]")
plt.legend()
plt.show()

