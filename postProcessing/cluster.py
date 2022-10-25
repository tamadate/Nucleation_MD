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
		self.con=con

	# function for the optimization
	def computeTth(self,t,v,x0):
		x=v*t+x0
		r=np.sum(x*x)
		fac=0
		if(t<0):
			fac=1e200
		return (self.con.delta2-r)**2+fac

	# get total calculation time
	def getTtot(self):
		self.U=np.loadtxt(self.con.directory+"U_"+str(self.con.I)+".dat") # Reading energy file(s)
		self.t_tot=np.max(self.U.T[0])*1e-15  # total simulation time [s]

	def Upot(self):
		self.getTtot()
		negs=np.where((self.U[:,0]>self.con.tEND*1e15))
		self.U=np.delete(self.U,negs,axis=0)
		self.plot.plotEnergies(self.U,self.con.teq,self.con.tEND,self.con.figOutput)

	def vapor(self):
		# Reading vapor in out time files
		inVapor=np.loadtxt(self.con.directory+"vapor_in_"+str(self.con.I)+".dat")
		negs=np.where((inVapor[:,1]>self.con.tEND*1e15))
		inVapor=np.delete(inVapor,negs,axis=0)
		Nin=np.size(inVapor.T[0])
		outVapor=np.loadtxt(self.con.directory+"vapor_out_"+str(self.con.I)+".dat")
		negs=np.where((outVapor[:,1]>self.con.tEND*1e15))
		outVapor=np.delete(outVapor,negs,axis=0)
		Nout=np.size(outVapor.T[0])

		# get total number of steps in analysis
		self.getTtot()
		Npost=int(self.t_tot/self.con.dt_post)

		times=np.arange(Npost)*self.con.dt_post
		Nstick=np.zeros(Npost)						# Number of vapors at the time
		usedLines=np.arange(Nout)					# Flag for outVapor (-1: used, 0: not used)
		ts=np.zeros(Nin)							# Residence time in MD
		tth=np.zeros(Nin)							# Residence time in theory

		for iin in np.arange(Nin):
			IDin=inVapor[iin][0]
			xin=inVapor[iin][5:8]
			vin=inVapor[iin][2:5]
			timeIn=inVapor[iin][1]*1e-15

			timeOut=self.t_tot
			# find same vapor ID
			for iout in np.arange(Nout):
				if(usedLines[iout]==-1):
					continue
				if(outVapor[iout][0]==IDin):
					timeOut=outVapor[iout][1]*1e-15
					usedLines[iout]=-1
					break
			Nstick[int(timeIn/self.con.dt_post):int(timeOut/self.con.dt_post)]+=1
			if(timeIn<self.con.teq):
				continue
			time3=minimize(self.computeTth,x0=100000,args=(xin,vin)).x[0]	# args=((vx,vy,vz),(x,y,z))
			ts[iin]=timeOut-timeIn-time3*1e-15	# time3 is theoretical residence time in interaction sphere
			tth[iin]=timeOut-timeIn	# time3 is theoretical residence time in interaction sphere

		self.plot.plotNvap(times,self.con.teq,Nstick,self.con.figOutput)
		self.plot.plotStickTimeDist(ts,tth,self.con.figOutput)

		#-----------------------------------------------------------------------------#
		# ignore if the regidence time is smaller than tcut
		negs=np.where((ts<self.con.tcut))
		tsCut=np.delete(ts,negs)
		tsave=np.average(tsCut)	# average sticking time [s]
		betaC=np.size(tsCut)/(self.t_tot-self.con.teq)	# vapor collision flux with ion [1/s]
		ram=betaC*tsave		# ramda for Poisson distribution

		NstickEq=Nstick[int(self.con.teq/self.con.dt_post):]
		Nbase=np.min(NstickEq)
		Nave=np.average(NstickEq)
		#ram=Nave-Nbase

		#ram-=Nbase-0.5
		nvap=np.arange(self.con.Nmax)
		ppoi=np.zeros(self.con.Nmax)
		psim=np.zeros(self.con.Nmax)
		for i in nvap:
			ppoi[i]=ram**i*np.exp(-ram)/math.factorial(i)
		for i in NstickEq:
		    psim[int(i)]+=1
		psim/=np.size(NstickEq)

		self.plot.plotStickVaporDist(nvap,ppoi,psim,Nbase,self.con.figOutput)
		f_sim=np.size(ts)/self.t_tot	# vapor flux into the interaction sphere [1/s]	# vapor flux into the interaction sphere in free molecular limit [1/s]
		print ("f_FM="+str(self.con.f_FM*1e-9)+"[1/ns]\tf_sim="+str(f_sim*1e-9)+"[1/ns]")


	def compute(self):
		self.Upot()
		self.vapor()
		self.plot.plotShow(self.con.figOutput,str(self.con.directory)+"fig"+str(self.con.I)+".png")
