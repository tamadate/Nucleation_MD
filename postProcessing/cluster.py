import numpy as np
import math
import os
import plot
import error
import conditions
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import statistics


class cluster:
	def __init__(self,con):
		os.system("make")

		self.plot=plot.plot()
		self.plot.pltNormal()
		self.plot.directory=con.directory
		self.con=con


	def Upot(self):
		## Reading energy file(s)
		Ubars,Ss,indexes,tENDs=[],[],[],[]
		for i in np.arange(self.con.pal):
			if(np.isin(i,self.con.error)):
				continue
			I=self.con.I+i
			U=np.loadtxt(self.con.directory+"U_"+str(I)+".dat")
			tENDs=np.append(tENDs,np.max(U.T[0]))  # total simulation time [s]
			negs=np.where((U[:,0]>self.con.tEND*1e15))
			U=np.delete(U,negs,axis=0)
			Ss=np.append(Ss,statistics.stdev(np.sum(U.T[1:],axis=0)))
			Ubars=np.append(Ubars,np.average(np.sum(U.T[1:],axis=0)))
			indexes=np.append(indexes,i)

		Umax=np.average(Ubars)+np.average(Ss)*2
		Umin=np.average(Ubars)-np.average(Ss)*2
		neg=np.where(((Ubars<Umin) | (Ubars>Umax)))

		if(np.size(neg)!=0):
			error.errorAdd(self.con,neg,indexes)
			Ubars,indexes,Ss,tENDs=np.delete(Ubars,neg),np.delete(indexes,neg),np.delete(Ss,neg),np.delete(tENDs,neg)
		plt.fill_between(indexes,Ubars-Ss*2,Ubars+Ss*2,color="gray",alpha=0.5,lw=0.01)
		plt.scatter(indexes,Ubars,c="black")
		plt.xlabel("Calculation id")
		plt.ylabel("Average energy [kcal/mol]")
		plt.show()
		plt.scatter(indexes,tENDs,c="black")
		plt.xlabel("Calculation id")
		plt.ylabel("Total calculation time [fs]")
		plt.show()

		self.plot.plotEnergies(U,self.con.teq,self.con.tEND,self.con.figOutput)

	def vapor(self):
		tss=np.zeros(0)
		psim=np.zeros(self.con.Nmax)
		for i in np.arange(self.con.pal):
			if(np.isin(i,self.con.error)):
				continue
			I=self.con.I+i

			## Reading vapor_in file
		    ## column[0]="vapor id" column[1]="entering time, tin"
			filePath=self.con.directory+"vapor_in_"+str(I)+".dat"
			if(os.stat(filePath).st_size == 0):	# if the file is empty
				continue
			inVapor=np.loadtxt(filePath)
			if(np.size(inVapor)<9):
				if(inVapor[1]>self.con.tEND*1e15): # ignore tin>tend
					continue
				inVapor=np.array([inVapor])
			else:
				negs=np.where((inVapor[:,1]>self.con.tEND*1e15))
				inVapor=np.delete(inVapor,negs,axis=0)

			## Reading vapor_out file
		    ## column[0]="vapor id" column[1]="leaving time, tin"
			outVapor=np.loadtxt(self.con.directory+"vapor_out_"+str(I)+".dat")
			if(np.size(outVapor)<9):
				if(outVapor[1]>self.con.tEND*1e15): # ignore tout>tend
					continue
				outVapor=np.array([outVapor])
			else:
				negs=np.where((outVapor.T[1]>self.con.tEND*1e15))
				outVapor=np.delete(outVapor,negs,axis=0)

			Npost=int(self.con.tEND/self.con.dt_post) # total number of steps in analysis
			times=np.arange(Npost)*self.con.dt_post
			Nstick=np.zeros(Npost)	# Number of vapors [N(t0),N(t0+self.con.dt_post),N(t1+2*self.con.dt_post),...,N(t_tot)]
			usedLines=np.arange(np.size(outVapor.T[0]))		# Flag for outVapor (equal to -1 if it is used)
			ts=np.zeros(np.size(inVapor.T[0]))	# Number of vapors [N(t0),N(t0+self.con.dt_post),N(t1+2*self.con.dt_post),...,N(t_tot)]
			tth=np.zeros(np.size(inVapor.T[0]))	# Number of vapors [N(t0),N(t0+self.con.dt_post),N(t1+2*self.con.dt_post),...,N(t_tot)]

			def function(t,v,x0):
				x=v*t+x0
				r=np.sum(x*x)
				fac=0
				if(t<0):
					fac=1e200
				return (self.con.delta2-r)**2+fac

			for iin in np.arange(np.size(inVapor.T[0])):
				time1=inVapor[iin][1]*1e-15
				time2=0
				for loop in np.arange(np.size(usedLines)):
					if(usedLines[loop]==-1):
						continue
					if(outVapor[loop][0]==inVapor[iin][0]):
						time2=outVapor[loop][1]*1e-15
						usedLines[loop]=-1
						break
				if(time2!=0):
					Nstick[int(time1/self.con.dt_post):int(time2/self.con.dt_post)]+=1
				else:
					Nstick[int(time1/self.con.dt_post):int(self.con.tEND/self.con.dt_post)]+=1
				result=minimize(function,x0=100000,args=(inVapor[iin][5:8],inVapor[iin][2:5]))	# args=((vx,vy,vz),(x,y,z))
				if(time1<self.con.teq):
					time1=self.con.tEND
				ts[iin]=time2-time1-result.x[0]*1e-15	# result.x[0] is theoretical residence time in interaction sphere
				tth[iin]=time2-time1	# result.x[0] is theoretical residence time in interaction sphere
			tss=np.append(tss,tth)
		#np.savetxt("ts.dat",ts)
		self.plot.plotNvap(times,self.con.teq,Nstick,self.con.figOutput)
		self.plot.plotStickTimeDist(tss,tth,self.con.figOutput)

		#-----------------------------------------------------------------------------#
		tsave=np.average(tss)	# average sticking time [s]
		betaC=np.size(np.delete(tss,negs))/(self.con.tEND-self.con.teq)/float(self.con.pal-np.size(self.con.error))	# vapor collision flux with ion [1/s]

		ram=betaC*tsave		# ramda for Poisson distribution

		Nbase=np.min(Nstick[int(self.con.teq/self.con.dt_post):])
		Nave=np.average(Nstick[int(self.con.teq/self.con.dt_post):])
		#ram=Nave-Nbase

		nv=np.arange(self.con.Nmax)
		ppoi=np.zeros(self.con.Nmax)

		for i in nv:
			ppoi[i]=ram**i*np.exp(-ram)/math.factorial(i)

		psim/=np.sum(psim)

		self.plot.plotStickVaporDist(nv,ppoi,psim,Nbase,self.con.figOutput)

		pv=self.con.pv0-Nbase 	# vapor pressure [Pa]
		C=pv/self.con.kb/self.con.T	# vapor concentraiton [1/m3]
		f_sim=np.size(ts)/self.con.tEND	# vapor flux into the interaction sphere [1/s]
		f_FM=self.con.delta2*np.pi*self.con.c*C	# vapor flux into the interaction sphere in free molecular limit [1/s]
		print ("f_FM="+str(f_FM*1e-9)+"[1/ns]\tf_sim="+str(f_sim*1e-9)+"[1/ns]")


	def compute(self):
		self.Upot()
		self.vapor()
