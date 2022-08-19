import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize

#-----------------------------------------------------------------------------#
#Cunningham
def calc_Cc(dp):
    ram=67e-9
    Kn=ram*2/dp
    A=1.165
    B=0.483
    C=0.997
    return 1+Kn*(A+B*np.exp(-C/Kn))

def pltNormal():
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'
    #plt.rcParams['figure.subplot.bottom'] = 0.2
    #plt.rcParams['figure.subplot.left'] = 0.2
    plt.rcParams['font.family'] = 'Arial'

def axNormal(ax):
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='x')
    ax.tick_params(axis='y')

pltNormal()
fig, axs = plt.subplots(2,2,figsize=(15,15))
for i in np.arange(4):
	axNormal(axs.flat[i])
#-----------------------------------------------------------------------------#


				## Reading energy file(s)
#-----------------------------------------------------------------------------#


U=np.loadtxt("../../../nucleation/Na/U_1.dat")
#K=np.loadtxt("K_1.dat")
labels=["Uion","Ugas","Uvap","Ugi","Ugg","Uvg","Uvi","Uvv"]

axs.flat[0].set_xlabel("Time [ns]")
axs.flat[0].set_ylabel("Energy [kcal/mol]")
for i in np.arange(np.size(U[0])-1):
	axs.flat[0].plot(U.T[0]*1e-6,U.T[i+1],label=labels[i])
axs.flat[0].plot(U.T[0]*1e-6,np.sum(U.T[1:],axis=0),label="Total")

#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
inVapor=np.loadtxt("../../../nucleation/Na/vapor_in_1.dat")
outVapor=np.loadtxt("../../../nucleation/Na/vapor_out_1.dat")

t_tot=0.8e-6	# total simulation time, t_tot [s]
dt=1e-15	# simulation time step, dt [s]
dt_post=1e-11								# dt in analysis, dt_post [s]
postStep=int(1e-11/dt)								# dt in analysis, dt_post [s]
Npost=int(t_tot/dt_post)					# total number of steps in analysis, Npost=/dt_post
teq=0.2e-6		# equilibliumed time

times=np.arange(Npost)*dt_post
Nstick=np.zeros(Npost)	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]
usedLines=np.arange(np.size(outVapor.T[0]))		# Flag for outVapor (equal to -1 if it is used)
ts=np.zeros(np.size(inVapor.T[0]))	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]

delta=1e-8	# diameter of interaction sphere [m]
M=0.018		# vapor mass [kg/mol]
kb=1.38e-23	# boltzmann constant [J/K]
R=8.314		# gas constant	[Jmol/K]
T=300.0		# temperature	[K]
c=(8*R*T/M/np.pi)**0.5	# vapor mean thermal speed [m/s]
pv=100		# vapor pressure [Pa]
C=pv/kb/T	# vapor concentraiton [1/m3]
f_FM=delta*delta*np.pi*c*C	# vapor flux into the interaction sphere in free molecular limit [1/s]

def function(t,v,x0):
	x=v*t+x0
	r=np.sum(x*x)
	fac=0
	if(t<0):
		fac=1e200
	return (delta*1e10-r)**2+fac

for iin in np.arange(np.size(inVapor.T[0])):
	time1=inVapor[iin][1]*dt
	if(time1<teq):
		time1=t_tot
	time2=0
	for loop in np.arange(np.size(usedLines)):
		if(usedLines[loop]==-1):
			continue
		if(outVapor[loop][0]==inVapor[iin][0]):
			time2=outVapor[loop][1]*dt
			usedLines[loop]=-1
			break

	if(time2!=0):
		Nstick[int(time1/dt_post):int(time2/dt_post)]+=1
	else:
		Nstick[int(time1/dt_post):int(t_tot/dt_post)]+=1
	result=minimize(function,x0=10000,args=(inVapor[iin][5:8],inVapor[iin][2:5]))	# args=((vx,vy,vz),(x,y,z))
	ts[iin]=time2-time1-result.x[0]*dt	# result.x[0] is theoretical residence time in interaction sphere

#np.savetxt("ts.dat",ts)

f_sim=np.size(ts)/t_tot	# vapor flux into the interaction sphere [1/s]
print ("f_FM="+str(f_FM*1e-9)+"[1/ns]\tf_sim="+str(f_sim*1e-9)+"[1/ns]")

negs=np.where(ts<1e-9)	# indexes of ts<1e-9
tsave=np.average(np.delete(ts,negs))	# average sticking time [s]
betaC=np.size(np.delete(ts,negs))/t_tot	# vapor collision flux with ion [1/s]
ram=betaC*tsave		# ramda for Poisson distribution

Nmax=10
nv=np.arange(Nmax)
ppoi=np.zeros(Nmax)
for i in np.arange(Nmax):
	ppoi[i]=ram**i*np.exp(-ram)/math.factorial(i)



axs.flat[1].set_xlabel("Time [ns]")
axs.flat[1].set_ylabel("Number of vapor in efective domain, Nvap [-]")
axs.flat[1].scatter(times*1e6,Nstick)

negs=np.where(ts<0)
axs.flat[2].set_xlabel("Logarithm of time [-]")
axs.flat[2].set_ylabel("Number of event [-]")
axs.flat[2].set_yscale("log")
axs.flat[2].hist(np.log10(np.delete(ts,negs)),bins=50)

#plt.savefig("figures/ONOFF_dist.png",dpi=1000)
plt.show()
