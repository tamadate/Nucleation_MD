import numpy as np
import math
import matplotlib.pyplot as plt
import os
import glob
import plot
from scipy.optimize import minimize

os.system("make")

#-----------------------------------------------------------------------------#
plot=plot.plot()

plot.pltNormal()
plot.fig, plot.axs = plt.subplots(1,2,figsize=(10,5))
for i in np.arange(2):
	plot.axNormal(plot.axs.flat[i])

plot.labelSize=15
plot.titleSize=18
plot.lineWidth=0.8
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

plot.plotMobilityShift(datas,directory)
