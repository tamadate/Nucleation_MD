import numpy as np
import math
import matplotlib.pyplot as plt
import os
import glob
from scipy.optimize import minimize

os.system("make")

                ##Cunningham
#-----------------------------------------------------------------------------#

def pltNormal():
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'
    #plt.rcParams['figure.subplot.bottom'] = 0.2
    #plt.rcParams['figure.subplot.left'] = 0.2
    #plt.rcParams['font.family'] = 'Arial'
    plt.rcParams["font.size"]=10

def axNormal(ax):
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='x')
    ax.tick_params(axis='y')

pltNormal()
fig, axs = plt.subplots(1,2,figsize=(10,5))
for i in np.arange(2):
	axNormal(axs.flat[i])

labelSize=15
titleSize=18
lineWidth=0.8

#-----------------------------------------------------------------------------#


				## Unknown parameters
#-----------------------------------------------------------------------------#
teq=0.2e-6     # equilibliumed time [s]
tcut=1e-13       # cut residence time [s]
Nmax=40         # Max number of sticking vapors
dt_post=1e-11   # dt in analysis, dt_post [s]
directory="../../../nucleation/NaCl/Na2Cl_toluene/"
#directory="../../../nucleation/angio2+/100/"
I=10

startTime=2e-9
endTime=5e-9

checkMode=0 # display
figOutput=1

#-----------------------------------------------------------------------------#
i=0
fileName=str(directory)+"DiffusionCoefficients."+str(i)
data=np.loadtxt(fileName,skiprows=1)
data=np.append([0],data)
datas=np.array([data])
for i in np.arange(100):
    fileName=str(directory)+"DiffusionCoefficients."+str(i)
    if(os.path.exists(fileName)):
        data=np.loadtxt(fileName,skiprows=1)
        if(math.isnan(data[0])):
            continue
        data=np.append([i],data)
        datas=np.append(datas,np.array([data]),axis=0)

print(datas)

#axs.flat[0].set_xlim([0,Xmax])
axs.flat[0].set_title("(a) Inverse of mobility shift",loc='left',fontsize=titleSize)
axs.flat[0].set_xlabel("Vapor pressure [Pa]",fontsize=labelSize)
axs.flat[0].set_ylabel(r"Mobility [cm$^2$/Vs]",fontsize=labelSize)
for data in datas:
    axs.flat[0].scatter(data[0]*10,(data[4]*data[6])**-1,color = "black")

axs.flat[1].set_title("(b) Mobility shift",loc='left',fontsize=titleSize)
axs.flat[1].set_xlabel("Vapor pressure [Pa]",fontsize=labelSize)
axs.flat[1].set_ylabel(r"Normalized mobility, $Z_p/Z_{p,0}$ [-]",fontsize=labelSize)
for data in datas:
    axs.flat[1].scatter(data[0]*10,(data[4]/datas[0][4]),color = "black")

plt.savefig(str(directory)+"diffusionSummary.png", dpi=1000)
plt.show()

#-----------------------------------------------------------------------------#
