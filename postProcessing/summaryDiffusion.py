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
plot.fig, plot.axs = plt.subplots(1,1,figsize=(5,5))
plot.axNormal(plot.axs)

plot.labelSize=15
plot.titleSize=18
plot.lineWidth=0.8
#-----------------------------------------------------------------------------#


				## Unknown parameters
#-----------------------------------------------------------------------------#
directory="/home/tama3rdgen/vaporUptake/bradykinin+2/v"
press=np.array((0,20,40,60,80,100,120,140,160,180,200))


#directory="/home/tama3rdgen/vaporUptake/angioII+2_new/v"
#press=np.array((0,20,60,80,100,120))


#-----------------------------------------------------------------------------#
i=0
datas=np.zeros((np.size(press),7))
for p in press:
	fileName=str(directory)+str(int(p/10))+"/"+str(int(p))+"Pa/DiffusionCoefficients.1"
	datas[i]=np.loadtxt(fileName,skiprows=1)
	i+=1
print(datas[0])
plot.plotMobilityShift(press,datas,directory)
