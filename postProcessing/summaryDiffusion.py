import numpy as np
import math
import matplotlib.pyplot as plt
import os
import glob
import plot
from scipy.optimize import minimize

os.system("make")

def pltNormal():
	plt.rcParams['ytick.direction'] = 'in'
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['figure.subplot.bottom'] = 0.2
	plt.rcParams['figure.subplot.left'] = 0.2
	#plt.rcParams['font.family'] = 'Arial'
	plt.rcParams["font.size"]=12

def axNormal(ax):
	ax.xaxis.set_ticks_position('both')
	ax.yaxis.set_ticks_position('both')
	ax.tick_params(axis='x')
	ax.tick_params(axis='y')

#-----------------------------------------------------------------------------#
pltNormal()

labelSize=15
titleSize=18
lineWidth=0.8

expLoad=np.loadtxt("experiment.csv",skiprows=1,delimiter=",")
#-----------------------------------------------------------------------------#
				## Unknown parameters
#-----------------------------------------------------------------------------#
directory="/media/tama3rdgen/6TB/vaporUptake/bradykinin/v"
press=np.array((0,20,40,60,80,100,120,140,180,240))
exp=np.array([[0,0,5,10,20,40,60,80,120,160,240,300],[1,1.0033,1.005275,1.001312,0.990309,0.987374,0.977268,0.979336,0.964985,0.96599,0.959359,0.951782]])
exp[1]*=0.59*2
z=2
molName="brad"
range=(0.55,0.60)


directory="/media/tama3rdgen/6TB/vaporUptake/angio2+_new/v"
press=np.array((0,20,40,60,80,100,120,140))
exp=np.array((expLoad.T[0],expLoad.T[11]))
z=2
molName="angio2"
range=(0.50,0.60)
'''
directory="/media/tama3rdgen/6TB/vaporUptake/angio1+_new/v"
press=np.array((0,10,20,40,60,80,140,180,240,300))
exp=np.array((expLoad.T[0],expLoad.T[9]))
z=1
molName="angio1"
range=(0.55,0.65)
'''
exp[0]*=111/300.0

#	Plot	
i=0
datas=np.zeros((np.size(press),7))
for p in press:
	fileName=str(directory)+str(int(p/10))+"/"+str(int(p))+"Pa/DiffusionCoefficients.1"
	datas[i]=np.loadtxt(fileName,skiprows=1)
	i+=1
fig, axs = plt.subplots(1,1,figsize=(5,5))
axNormal(axs)

axs.set_xlabel("Vapor pressure [Pa]",fontsize=labelSize)
axs.set_ylabel(r"Reduced mobility, $Z_p/z$ [cm$^2$ V$^{-1}$ s$^{-1}$]",fontsize=labelSize)
coeff=1.6e-19/1.38e-23/300.0
axs.scatter(press,datas.T[3]*coeff,color = "black",label="MSD")
axs.scatter(press,datas.T[4]*coeff,color = "red",label="VACF")
axs.scatter(exp[0],exp[1]/float(z),facecolor="none",edgecolor="blue",marker="^",label="Experiment")
axs.set_xlim(-10,150)
axs.set_ylim(range)
axs.legend()
plt.savefig(str(directory)+"diffusionSumABS_"+molName+".png", dpi=1000)

fig, axs = plt.subplots(1,1,figsize=(5,5))
axNormal(axs)
exp[1]/=exp[1][0]

axs.set_xlabel("Vapor pressure [Pa]",fontsize=labelSize)
axs.set_ylabel(r"Mobility shift, $Z_p/Z_{p,0}$ [-]",fontsize=labelSize)
K0=(datas.T[3][0]+datas.T[4][0])*0.5
axs.scatter(press,datas.T[3]/K0,color = "black",label="MSD")
axs.scatter(press,datas.T[4]/K0,color = "red",label="VACF")
axs.scatter(exp[0],exp[1],facecolor="none",edgecolor="blue",marker="^",label="Experiment")
axs.set_xlim(-10,150)
axs.set_ylim(0.9,1.05)
axs.legend()
plt.savefig(str(directory)+"diffSum_"+molName+".png", dpi=1000)

plt.show()

