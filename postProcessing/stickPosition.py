import numpy as np
import math
import os
import plot
import conditions
import matplotlib.pyplot as plt
from scipy.optimize import minimize


class stickPosition:
	def __init__(self,con):
		os.system("make")

		self.plot=plot.plot()
		self.plot.pltNormal()
		self.plot.fig, self.plot.axs = plt.subplots(1,2,figsize=(12,5))
		for i in np.arange(2):
			self.plot.axNormal(self.plot.axs.flat[i])
		self.con=con

	def compute(self):
		amino=np.loadtxt(self.con.inputDirectory+"amino.loc",dtype=np.str)
		vaporID=np.loadtxt(self.con.inputDirectory+"vaporID.loc",dtype=np.int)-1
		aminoNames=np.unique(amino)
		dist=np.zeros(np.size(vaporID))
		distAmino=np.zeros(np.size(aminoNames))
		label=[]
		Ntot=0.0
		for i in np.arange(np.size(vaporID)):
			label.append(str(int(vaporID[i])+1))
		for rep in self.con.replica:
			#data=np.loadtxt(self.con.directory+str(rep)+"/stickPositionLog_"+str(self.con.I)+".dat")
			data=np.loadtxt(self.con.directory+str(rep)+"/stickPositionLog_"+str(rep)+".dat")
			negs=np.where(data.T[3]>5)
			data=np.delete(data,negs,axis=0)
			for i in np.arange(np.size(vaporID)):
				gets=np.where(data.T[2]==vaporID[i])
				N=np.size(data[gets].T[0])
				dist[i]+=N
				distAmino[np.where(aminoNames==amino[i])]+=N
				Ntot+=N
		self.plot.plotStickLocation(label,dist/Ntot*100,aminoNames,distAmino/Ntot*100,self.con.directory)
