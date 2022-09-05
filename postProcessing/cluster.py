import numpy as np
import math
import matplotlib.pyplot as plt
import os
from scipy.optimize import minimize

os.system("make")

                ##Cunningham
#-----------------------------------------------------------------------------#

def pltNormal():
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'
    #plt.rcParams['figure.subplot.bottom'] = 0.2
    #plt.rcParams['figure.subplot.left'] = 0.2
    #plt.rcParams['font.family'] = 'Arial'
    plt.rcParams["font.size"]=8

def axNormal(ax):
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='x')
    ax.tick_params(axis='y')

pltNormal()
fig, axs = plt.subplots(3,2,figsize=(8,10))
for i in np.arange(6):
	axNormal(axs.flat[i])

labelSize=10
titleSize=12
lineWidth=0.8

#-----------------------------------------------------------------------------#


				## Unknown parameters
#-----------------------------------------------------------------------------#
teq=1e-6     # equilibliumed time [s]
tcut=1e-13       # cut residence time [s]
Nmax=40         # Max number of sticking vapors
dt_post=1e-11   # dt in analysis, dt_post [s]
directory="../../../nucleation/NaCl/Na2Cl/"
#directory="../../../nucleation/angio2+/100/"
I=3

startTime=2e-9
endTime=5e-9

checkMode=0
figOutput=1

#-----------------------------------------------------------------------------#



				## Reading energy file(s)
#-----------------------------------------------------------------------------#

U=np.loadtxt(str(directory)+"U_"+str(I)+".dat")
#K=np.loadtxt("K_1.dat")
labels=["Uion","Ugas","Uvap","Ugi","Ugg","Uvg","Uvi","Uvv"]
display=[1,0,0,0,0,0,1,1]
colors=["blue","blue","blue","blue","blue","blue","aqua","navy"]

t_tot=np.max(U.T[0])*1e-15  # total simulation time [s]
Xmax=int(t_tot/200e-9)*200
Xmax*=1.4

#axs.flat[0].set_xlim([0,Xmax])
axs.flat[0].set_title("(a) Potential Energies",loc='left',fontsize=titleSize)
axs.flat[0].set_xlabel("Time [ns]",fontsize=labelSize)
axs.flat[0].set_ylabel("Energy [kcal/mol]",fontsize=labelSize)
for i in np.arange(np.size(U[0])-1):
    if(display[i]):
	       axs.flat[0].plot(U.T[0]*1e-6,U.T[i+1],label=labels[i],linewidth=lineWidth,color=colors[i])
axs.flat[0].plot(U.T[0]*1e-6,np.sum(U.T[1:],axis=0),label="Total",linewidth=lineWidth,color="black")
axs.flat[0].axvline(x = teq*1e9, ls='--', color = 'black',linewidth=lineWidth)
axs.flat[0].legend()

#-----------------------------------------------------------------------------#



				## Reading vapor in out time files
#-----------------------------------------------------------------------------#

inVapor=np.loadtxt(str(directory)+"vapor_in_"+str(I)+".dat")
outVapor=np.loadtxt(str(directory)+"vapor_out_"+str(I)+".dat")

Npost=int(t_tot/dt_post)					# total number of steps in analysis, Npost=/dt_post

times=np.arange(Npost)*dt_post
Nstick=np.zeros(Npost)	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]
usedLines=np.arange(np.size(outVapor.T[0]))		# Flag for outVapor (equal to -1 if it is used)
ts=np.zeros(np.size(inVapor.T[0]))	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]
tth=np.zeros(np.size(inVapor.T[0]))	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]

delta=1e-8	# diameter of interaction sphere [m]
delta2=(delta*1e10)**2      # square of delta [ang^2]
M=1/(1/0.032)		# vapor mass [kg/mol]
kb=1.38e-23	# boltzmann constant [J/K]
R=8.314		# gas constant	[Jmol/K]
T=300.0		# temperature	[K]
c=(8*R*T/M/np.pi)**0.5	# vapor mean thermal speed [m/s]

def function(t,v,x0):
	x=v*t+x0
	r=np.sum(x*x)
	fac=0
	if(t<0):
		fac=1e200
	return (delta2-r)**2+fac

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
        Nstick[int(time1/dt_post):int(time2/dt_post)]+=1
    else:
        Nstick[int(time1/dt_post):int(t_tot/dt_post)]+=1
    result=minimize(function,x0=100000,args=(inVapor[iin][5:8],inVapor[iin][2:5]))	# args=((vx,vy,vz),(x,y,z))
    if(time1<teq):
        time1=t_tot
    ts[iin]=time2-time1-result.x[0]*1e-15	# result.x[0] is theoretical residence time in interaction sphere
    tth[iin]=time2-time1	# result.x[0] is theoretical residence time in interaction sphere

#np.savetxt("ts.dat",ts)

f_sim=np.size(ts)/t_tot	# vapor flux into the interaction sphere [1/s]

axs.flat[1].set_title("(b) MeOHs in effective domain",loc='left',fontsize=titleSize)
axs.flat[1].set_xlabel("Time [ns]",fontsize=labelSize)
axs.flat[1].set_ylabel("$\it {N}$$_ {vap}$ [-]",fontsize=labelSize)
axs.flat[1].axvline(x = teq*1e9, color = 'black', ls="--",linewidth=lineWidth)
axs.flat[1].scatter(times*1e9,Nstick,color="black",s=20)

negs=np.where(ts<=0)
axs.flat[2].set_title("(c) Vapor sticking time distribution",loc='left',fontsize=titleSize)
axs.flat[2].set_xlabel("Logarithm of time [-]",fontsize=labelSize)
axs.flat[2].set_ylabel("Number of event [-]",fontsize=labelSize)
axs.flat[2].set_yscale("log")
axs.flat[2].hist(np.log10(np.delete(ts,negs)),alpha=0.3,bins=50,label="$\it {t}$$_{sim}$",color="blue")
axs.flat[2].hist(np.log10(np.delete(tth,negs)),alpha=0.3,bins=30,label="$\it {t}$$_ {sim}$-$\it t$$_ {th}$",color="black")
#axs.flat[2].axvline(x = np.log10(tcut), color = 'black', ls="--",linewidth=lineWidth)
axs.flat[2].legend()

if(checkMode==1):
    plt.show()

#-----------------------------------------------------------------------------#




#-----------------------------------------------------------------------------#
negs=np.where((ts<tcut))	# indexes of ts<tcut
tsave=np.average(np.delete(ts,negs))	# average sticking time [s]
betaC=np.size(np.delete(ts,negs))/(t_tot-teq)	# vapor collision flux with ion [1/s]
ram=betaC*tsave		# ramda for Poisson distribution

Nbase=np.min(Nstick[int(teq/dt_post):])
Nave=np.average(Nstick[int(teq/dt_post):])
#ram=Nave-Nbase

nv=np.arange(Nmax)
ppoi=np.zeros(Nmax)
psim=np.zeros(Nmax)
for i in nv:
	ppoi[i]=ram**i*np.exp(-ram)/math.factorial(i)
for i in Nstick[int(teq/dt_post):]:
    psim[int(i)]+=1
psim/=np.size(Nstick[int(teq/dt_post):])


axs.flat[3].set_title("(d) Sticking vapor distribution",loc='left',fontsize=titleSize)
axs.flat[3].set_title("")
axs.flat[3].set_xlabel("Number of sticking vapor",fontsize=labelSize)
axs.flat[3].set_ylabel("Frequency [-]",fontsize=labelSize)
axs.flat[3].set_yscale("log")
axs.flat[3].set_ylim([1e-5,1])
axs.flat[3].scatter(nv,ppoi,label="Poisson",s=20,color="black")
axs.flat[3].scatter(nv,psim,label="Simulation",s=20,marker="v",color="blue")
axs.flat[3].scatter(nv-Nbase,psim,label="Simulation("+str(int(-Nbase))+")",s=20,marker="^",color="cyan")
axs.flat[3].legend()

pv=10*I-Nbase 	# vapor pressure [Pa]
C=pv/kb/T	# vapor concentraiton [1/m3]
f_FM=delta*delta*np.pi*c*C	# vapor flux into the interaction sphere in free molecular limit [1/s]
print ("f_FM="+str(f_FM*1e-9)+"[1/ns]\tf_sim="+str(f_sim*1e-9)+"[1/ns]")

#-----------------------------------------------------------------------------#


ionData=np.loadtxt(str(directory)+"ion_300_"+str(I)+".dat")
dtInput=ionData[1][0]-ionData[0][0]
os.system("./mobility.out "+str(I)+" 1 "+str(teq)+" "+str(startTime)+" "+str(endTime)+" "+str(dtInput*1e-9)+" "+str(directory))
#"ion_300_%d.dat"

MSDVAF=np.loadtxt(str(directory)+"TIME_MSD_VAF."+str(I))
diffusionData=np.loadtxt(str(directory)+"DiffusionCoefficients."+str(I),skiprows=1)
axs.flat[4].set_title("(e) Normalized MSD",loc='left',fontsize=titleSize)
axs.flat[4].set_xlabel("Time [ns]",fontsize=labelSize)
axs.flat[4].set_ylabel("Normalized MSD [-]",fontsize=labelSize)
axs.flat[4].scatter(MSDVAF.T[0],MSDVAF.T[1]/np.max(MSDVAF.T[1]),color="black",s=10)
axs.flat[4].text(1, 0.8, "K = "+'{:.2f}'.format(diffusionData[3]*diffusionData[5])+" cm$^ 2$/Vs", fontsize = 10)

axs.flat[5].set_title("(f) Normalized VAF",loc='left',fontsize=titleSize)
axs.flat[5].set_xlabel("Time [ns]",fontsize=labelSize)
axs.flat[5].set_ylabel("Normalized VAF [-]",fontsize=labelSize)
axs.flat[5].scatter(MSDVAF.T[0],MSDVAF.T[2]/MSDVAF.T[2][0],color="black",s=10)
axs.flat[5].text(2, 0.8, "K = "+'{:.2f}'.format(diffusionData[4]*diffusionData[5])+" cm$^ 2$/Vs", fontsize = 10)

#-----------------------------------------------------------------------------#

fig.tight_layout()

if(figOutput):
    plt.savefig(str(directory)+"fig"+str(I)+".png", dpi=1000)

plt.show()
