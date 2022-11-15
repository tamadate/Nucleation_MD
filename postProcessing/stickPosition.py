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
		amino=np.loadtxt(self.con.directory+"input/amino.loc",dtype=np.str)
		vaporIDstr=np.loadtxt(self.con.directory+"input/vaporID.loc",dtype=np.str)
		vaporID=np.loadtxt(self.con.directory+"input/vaporID.loc")
		aminoNames=np.unique(amino)
		dist=np.zeros(np.size(vaporID))
		distAmino=np.zeros(np.size(aminoNames))
		for j in np.arange(self.con.pal):
			if(np.isin(j,self.con.error)):
				continue
			I=self.con.I+j
			file=self.con.directory+"stickPositionLog_"+str(I)+".dat"
			if(os.path.exists(file)==False):
				continue
			data=np.loadtxt(file)
			#negs=np.where(data.T[3]>20)
			#data=np.delete(data,negs,axis=0)
			for i in np.arange(np.size(vaporID)):
				gets=np.where(data.T[2]+1==vaporID[i])
				N=np.size(data[gets].T[0])
				dist[i]+=N
				distAmino[np.where(aminoNames==amino[i])]+=N

		self.plot.plotStickLocation(vaporIDstr,dist,aminoNames,distAmino,self.con.directory)
