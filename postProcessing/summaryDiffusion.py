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
directory="/media/tama3rdgen/6TB/vaporUptake/bradykinin/v"
press=np.array((0,20,40,60,80,100,120,140,180,240))

#directory="/media/tama3rdgen/6TB/vaporUptake/angio2+_new/v"
#press=np.array((0,20,40,60,80,100,120,140))

#directory="/media/tama3rdgen/6TB/vaporUptake/angio1+_new/v"
#press=np.array((0,10,20,40,60,80,140,180,240,300))

exp=np.array([[0,0,5,10,20,40,60,80,120,160,240,300],[1,1.0033,1.005275,1.001312,0.990309,0.987374,0.977268,0.979336,0.964985,0.96599,0.959359,0.951782]])
exp[0]*=111/300.0
exp[1]=exp[1]*0+0.6

#exp=np.array([[0],[1]])


#-----------------------------------------------------------------------------#
i=0
datas=np.zeros((np.size(press),7))
for p in press:
	fileName=str(directory)+str(int(p/10))+"/"+str(int(p))+"Pa/DiffusionCoefficients.1"
	datas[i]=np.loadtxt(fileName,skiprows=1)
	i+=1
print(datas[0])
datas[0]*=1.0
plot.plotMobilityShift(press,datas,directory,exp)
