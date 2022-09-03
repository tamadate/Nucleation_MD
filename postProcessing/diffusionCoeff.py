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
fig, axs = plt.subplots(1,2,figsize=(15,5))
for i in np.arange(2):
	axNormal(axs.flat[i])
#-----------------------------------------------------------------------------#


				## Unknown parameters
#-----------------------------------------------------------------------------#
teq=0.2e-6     # equilibliumed time [s]
tcut=1e-13       # cut residence time [s]
Nmax=30         # Max number of sticking vapors
dt=1e-15	    # simulation time step, dt [s]
dt_post=1e-11   # dt in analysis, dt_post [s]
directory="../../../nucleation/valine/"
#directory="./"
I=10
checkMode=0

startTime=2e-9
endTime=5e-9
#-----------------------------------------------------------------------------#



ionData=np.loadtxt(str(directory)+"ion_300_"+str(I)+".dat")
dtInput=ionData[1][0]-ionData[0][0]
os.system("./mobility.out "+str(I)+" 1 "+str(teq)+" "+str(startTime)+" "+str(endTime)+" "+str(dtInput*1e-9)+" "+str(directory)+"ion_300_%d.dat")
#"ion_300_%d.dat"

MSDVAF=np.loadtxt("TIME_MSD_VAF."+str(I))
axs.flat[0].set_xlabel("Time [ns]",fontsize=20)
axs.flat[0].set_ylabel("MSD [ang$^ 2$]",fontsize=20)
axs.flat[0].scatter(MSDVAF.T[0],MSDVAF.T[1]/np.max(MSDVAF.T[1]))

axs.flat[1].set_xlabel("Time [ns]",fontsize=20)
axs.flat[1].set_ylabel("VAF [ang$^ 2$]",fontsize=20)
axs.flat[1].scatter(MSDVAF.T[0],MSDVAF.T[2]/MSDVAF.T[2][0])


plt.show()
