import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize

                ##Cunningham
#-----------------------------------------------------------------------------#

def pltNormal():
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'
    #plt.rcParams['figure.subplot.bottom'] = 0.2
    #plt.rcParams['figure.subplot.left'] = 0.2
    #plt.rcParams['font.family'] = 'Arial'
    plt.rcParams["font.size"]=12

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


				## Unknown parameters
#-----------------------------------------------------------------------------#
teq=0.1e-6     # equilibliumed time [s]
tcut=1e-9       # cut residence time [s]
Nmax=30         # Max number of sticking vapors
dt=1e-15	    # simulation time step, dt [s]
dt_post=1e-11   # dt in analysis, dt_post [s]
directory="../../../nucleation/lysine/"
#directory="../../../nucleation/angio2+/100/"
I=10
checkMode=0
#-----------------------------------------------------------------------------#



				## Reading energy file(s)
#-----------------------------------------------------------------------------#

U=np.loadtxt(str(directory)+"U_"+str(I)+".dat")
#K=np.loadtxt("K_1.dat")
labels=["Uion","Ugas","Uvap","Ugi","Ugg","Uvg","Uvi","Uvv"]

t_tot=np.max(U.T[0])*1e-15  # total simulation time [s]
Xmax=int(t_tot/200e-9)*200
Xmax*=1.4

#axs.flat[0].set_xlim([0,Xmax])
axs.flat[0].set_xlabel("Time [ns]",fontsize=20)
axs.flat[0].set_ylabel("Energy [kcal/mol]",fontsize=20)
for i in np.arange(np.size(U[0])-1):
	axs.flat[0].plot(U.T[0]*1e-6,U.T[i+1],label=labels[i])
axs.flat[0].plot(U.T[0]*1e-6,np.sum(U.T[1:],axis=0),label="Total")
axs.flat[0].axvline(x = teq*1e9, ls='--', color = 'black')
axs.flat[0].axvline(x = t_tot*1e9, ls='--', color = 'black')
axs.flat[0].legend()

#-----------------------------------------------------------------------------#



				## Reading vapor in out time files
#-----------------------------------------------------------------------------#

inVapor=np.loadtxt(str(directory)+"vapor_in_"+str(I)+".dat")
outVapor=np.loadtxt(str(directory)+"vapor_out_"+str(I)+".dat")

postStep=int(dt_post/dt)								# dt in analysis, dt_post [s]
Npost=int(t_tot/dt_post)					# total number of steps in analysis, Npost=/dt_post

times=np.arange(Npost)*dt_post
Nstick=np.zeros(Npost)	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]
usedLines=np.arange(np.size(outVapor.T[0]))		# Flag for outVapor (equal to -1 if it is used)
ts=np.zeros(np.size(inVapor.T[0]))	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]
tth=np.zeros(np.size(inVapor.T[0]))	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]

delta=1e-8	# diameter of interaction sphere [m]
delta2=(delta*1e10)**2      # square of delta [ang^2]
M=1/(1/0.018+1/0.023)		# vapor mass [kg/mol]
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
	return (delta2-r)**2+fac

for iin in np.arange(np.size(inVapor.T[0])):
    time1=inVapor[iin][1]*dt
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
    result=minimize(function,x0=100000,args=(inVapor[iin][5:8],inVapor[iin][2:5]))	# args=((vx,vy,vz),(x,y,z))
    if(time1<teq):
        time1=t_tot
    ts[iin]=time2-time1-result.x[0]*dt	# result.x[0] is theoretical residence time in interaction sphere
    tth[iin]=time2-time1	# result.x[0] is theoretical residence time in interaction sphere

#np.savetxt("ts.dat",ts)

f_sim=np.size(ts)/t_tot	# vapor flux into the interaction sphere [1/s]
print ("f_FM="+str(f_FM*1e-9)+"[1/ns]\tf_sim="+str(f_sim*1e-9)+"[1/ns]")


axs.flat[1].set_xlabel("Time [ns]",fontsize=20)
axs.flat[1].set_ylabel("$\it {N}$$_ {vap}$ [-]",fontsize=20)
axs.flat[1].axvline(x = teq*1e9, color = 'black', ls="--")
axs.flat[1].scatter(times*1e9,Nstick)

negs=np.where(ts<=0)
axs.flat[2].set_xlabel("Logarithm of time [-]",fontsize=20)
axs.flat[2].set_ylabel("Number of event [-]",fontsize=20)
axs.flat[2].set_yscale("log")
axs.flat[2].hist(np.log10(np.delete(ts,negs)),alpha=0.3,bins=50,label="$\it {t}$$_{sim}$")
axs.flat[2].hist(np.log10(np.delete(tth,negs)),alpha=0.3,bins=30,label="$\it {t}$$_ {sim}$-$\it t$$_ {th}$")
axs.flat[2].axvline(x = np.log10(tcut), color = 'black', ls="--")
axs.flat[2].legend()

if(checkMode==1):
    plt.show()

#-----------------------------------------------------------------------------#




#-----------------------------------------------------------------------------#
negs=np.where(ts<tcut)	# indexes of ts<1e-9
tsave=np.average(np.delete(ts,negs))	# average sticking time [s]
betaC=np.size(np.delete(ts,negs))/t_tot	# vapor collision flux with ion [1/s]
ram=betaC*tsave		# ramda for Poisson distribution

nv=np.arange(Nmax)
ppoi=np.zeros(Nmax)
psim=np.zeros(Nmax)
for i in nv:
	ppoi[i]=ram**i*np.exp(-ram)/math.factorial(i)
for i in Nstick[int(teq/dt_post):]:
    psim[int(i)]+=1
psim/=np.size(Nstick[int(teq/dt_post):])


axs.flat[3].set_xlabel("Number of sticking vapor",fontsize=20)
axs.flat[3].set_ylabel("Frequency [-]",fontsize=20)
axs.flat[3].set_yscale("log")
axs.flat[3].set_ylim([1e-5,1])
axs.flat[3].scatter(nv,ppoi,label="Poisson")
axs.flat[3].scatter(nv,psim,label="Simulation")
axs.flat[3].scatter(nv-8,psim,label="Simulation(-8)")
axs.flat[3].legend()

#-----------------------------------------------------------------------------#

plt.savefig("fig.png", dpi=1000)
plt.show()
