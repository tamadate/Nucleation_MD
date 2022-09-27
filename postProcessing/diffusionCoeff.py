import numpy as np
import math
import matplotlib.pyplot as plt
import os
from scipy.optimize import minimize

                ##Cunningham
#-----------------------------------------------------------------------------#

def pltNormal():
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'
    #plt.rcParams['figure.subplot.bottom'] = 0.2
    #plt.rcParams['figure.subplot.left'] = 0.2
    #plt.rcParams['font.family'] = 'Arial'
    plt.rcParams["font.size"]=12

def axNormal(ax):
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='x')
    ax.tick_params(axis='y')

pltNormal()
fig, axs = plt.subplots(2,2,figsize=(15,10))
for i in np.arange(4):
	axNormal(axs.flat[i])
#-----------------------------------------------------------------------------#


				## Unknown parameters
#-----------------------------------------------------------------------------#
teq=0.01e-6     # equilibliumed time [s]
tEND=0.8e-6
Nmax=30         # Max number of sticking vapors
dt=0.5e-15	    # simulation time step, dt [s]
directory="/home/tama3rdgen/vaporUptake/angioII+1/withoutVap/"
I=0

startTime=2e-9
endTime=5e-9
#-----------------------------------------------------------------------------#

os.system("make")

ionData=np.loadtxt(str(directory)+"ion_300_"+str(I)+".dat")
dtInput=ionData[1][0]-ionData[0][0]
os.system("./mobility.out "+str(I)+" 1 "+str(teq)+" "+str(startTime)+" "+str(endTime)+" "+str(dtInput*1e-9)+" "+str(directory)+" "+str(tEND))
#"ion_300_%d.dat"

MSDVAF=np.loadtxt(str(directory)+"TIME_MSD_VAF."+str(I))
axs.flat[0].set_xlabel("Time [ns]",fontsize=20)
axs.flat[0].set_ylabel("MSD [ang$^ 2$]",fontsize=20)
axs.flat[0].scatter(MSDVAF.T[0],MSDVAF.T[1])

axs.flat[1].set_xlabel("Time [ns]",fontsize=20)
axs.flat[1].set_ylabel("VAF [ang$^ 2$]",fontsize=20)
axs.flat[1].scatter(MSDVAF.T[0],MSDVAF.T[2])


plt.show()
