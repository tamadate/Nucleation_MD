import numpy as np
import math
import matplotlib.pylab as plt
from scipy.optimize import minimize

inVapor=np.loadtxt("vapor_in_1.dat")
outVapor=np.loadtxt("vapor_out_1.dat")

t_tot=0.8e-6	# total simulation time, t_tot [s]
dt=1e-15	# simulation time step, dt [s]
dt_post=1e-11								# dt in analysis, dt_post [s]
postStep=int(1e-11/dt)								# dt in analysis, dt_post [s]
Npost=int(t_tot/dt_post)					# total number of steps in analysis, Npost=/dt_post
teq=0.2e-6

times=np.arange(Npost)*dt_post
Nstick=np.zeros(Npost)	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]
usedLines=np.arange(np.size(outVapor.T[0]))		# Flag for outVapor (equal to -1 if it is used)
ts=np.zeros(np.size(inVapor.T[0]))	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]

def function(t,v,x0):
	r0=100
	x=v*t+x0
	r=np.sum(x*x)
	fac=0
	if(t<0):
		fac=1e200
	return (r0-r)**2+fac

iin=0
for i in inVapor:
	time1=i[1]*dt
	if(time1<teq):
		time1=t_tot
	time2=0
	for loop in np.arange(np.size(usedLines)):
		if(usedLines[loop]==-1):
			continue
		if(outVapor[loop][0]==i[0]):
			time2=outVapor[loop][1]*dt
			usedLines[loop]=-1
			break
	
	result=minimize(function,x0=10000,args=(i[5:8],i[2:5]))
	ts[iin]=time2-time1#-result.x[0]*dt
	iin+=1
	


np.savetxt("ts.dat",ts)

negs=np.where(ts<1e-9)
tsave=np.average(np.delete(ts,negs))
betaC=np.size(np.delete(ts,negs))/t_tot

M=0.018
kb=1.38e-23
R=8.314
T=300.0
delta=1e-8
c=(8*R*T/M/np.pi)**0.5
print (c)
beta=delta*delta*np.pi*c

pv=100
C=pv/kb/T

ram=betaC*tsave
print(ram)

nv=np.arange(10)
ppoi=np.zeros(10)
for i in np.arange(10):
	ppoi[i]=ram**i*np.exp(-ram)/math.factorial(i)

print(ppoi)

negs=np.where(ts<0)
plt.hist(np.log10(np.delete(ts,negs)),bins=50)
plt.yscale("log")
plt.xlabel("Time [ns]")
plt.ylabel("Number of vapor in efective domain, Nvap [-]")
plt.show()

