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
		data=np.loadtxt(self.con.directory+"stickPositionLog_"+str(self.con.I)+".dat")
		amino=np.loadtxt(self.con.directory+"input/amino.loc",dtype=np.str)
		negs=np.where(data.T[3]>20)
		data=np.delete(data,negs,axis=0)
		vaporID=np.unique(data.T[2])
		aminoNames=np.unique(amino)
		dist=np.zeros(np.size(vaporID))
		distAmino=np.zeros(np.size(aminoNames))
		label=[]
		for i in np.arange(np.size(vaporID)):
			gets=np.where(data.T[2]==vaporID[i])
			N=np.size(data[gets].T[0])
			dist[i]+=N
			distAmino[np.where(aminoNames==amino[i])]+=N
			label.append(str(int(vaporID[i])))
		self.plot.plotStickLocation(label,dist/np.size(data.T[0])*100,aminoNames,distAmino/np.size(data.T[0])*100,self.con.directory)
