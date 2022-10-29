import numpy as np
import math
import os
import plot
import conditions
import matplotlib.pyplot as plt
from scipy.optimize import minimize


class cluster:
	def __init__(self,con):
		os.system("make")

		self.plot=plot.plot()
		self.plot.pltNormal()
		self.plot.directory=con.directory
		self.con=con


	def Upot(self):
		## Reading energy file(s)
		U=np.loadtxt(self.con.directory+"U_"+str(self.con.I)+".dat")
		self.t_tot=np.max(U.T[0])*1e-15  # total simulation time [s]
		negs=np.where((U[:,0]>self.con.tEND*1e15))
		U=np.delete(U,negs,axis=0)
		self.plot.plotEnergies(U,self.con.teq,self.con.tEND,0)


	def vapor(self):
		## Reading vapor in out time files
		inVapor=np.loadtxt(self.con.directory+"vapor_in_"+str(self.con.I)+".dat")
		negs=np.where((inVapor[:,1]>self.con.tEND*1e15))
		inVapor=np.delete(inVapor,negs,axis=0)
		outVapor=np.loadtxt(self.con.directory+"vapor_out_"+str(self.con.I)+".dat")
		negs=np.where((outVapor[:,1]>self.con.tEND*1e15))
		outVapor=np.delete(outVapor,negs,axis=0)
		Npost=int(self.t_tot/self.con.dt_post)					# total number of steps in analysis, Npost=/self.con.dt_post

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
		        Nstick[int(time1/self.con.dt_post):int(self.t_tot/self.con.dt_post)]+=1
		    result=minimize(function,x0=100000,args=(inVapor[iin][5:8],inVapor[iin][2:5]))	# args=((vx,vy,vz),(x,y,z))
		    if(time1<self.con.teq):
		        time1=self.t_tot
		    ts[iin]=time2-time1-result.x[0]*1e-15	# result.x[0] is theoretical residence time in interaction sphere
		    tth[iin]=time2-time1	# result.x[0] is theoretical residence time in interaction sphere

		#np.savetxt("ts.dat",ts)

		self.plot.plotNvap(times,self.con.teq,Nstick,self.con.figOutput)
		self.plot.plotStickTimeDist(ts,tth,self.con.figOutput)



		#-----------------------------------------------------------------------------#
		negs=np.where((ts<self.con.tcut))	# indexes of ts<self.con.tcut
		tsave=np.average(np.delete(ts,negs))	# average sticking time [s]
		betaC=np.size(np.delete(ts,negs))/(self.t_tot-self.con.teq)	# vapor collision flux with ion [1/s]
		ram=betaC*tsave		# ramda for Poisson distribution

		Nbase=np.min(Nstick[int(self.con.teq/self.con.dt_post):])
		Nave=np.average(Nstick[int(self.con.teq/self.con.dt_post):])
		#ram=Nave-Nbase

		nv=np.arange(self.con.Nmax)
		ppoi=np.zeros(self.con.Nmax)
		psim=np.zeros(self.con.Nmax)
		for i in nv:
			ppoi[i]=ram**i*np.exp(-ram)/math.factorial(i)
		for i in Nstick[int(self.con.teq/self.con.dt_post):]:
		    psim[int(i)]+=1
		psim/=np.size(Nstick[int(self.con.teq/self.con.dt_post):])

		self.plot.plotStickVaporDist(nv,ppoi,psim,Nbase,self.con.figOutput)


		pv=self.con.pv0-Nbase 	# vapor pressure [Pa]
		C=pv/self.con.kb/self.con.T	# vapor concentraiton [1/m3]
		f_sim=np.size(ts)/self.t_tot	# vapor flux into the interaction sphere [1/s]
		f_FM=self.con.delta2*np.pi*self.con.c*C	# vapor flux into the interaction sphere in free molecular limit [1/s]
		print ("f_FM="+str(f_FM*1e-9)+"[1/ns]\tf_sim="+str(f_sim*1e-9)+"[1/ns]")


	def compute(self):
		self.Upot()
		self.vapor()
