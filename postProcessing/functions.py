import numpy as np
import os
import plot
from scipy.optimize import minimize



pv=pv0-Nbase 	# vapor pressure [Pa]
C=pv/kb/T	# vapor concentraiton [1/m3]
f_sim=np.size(ts)/t_tot	# vapor flux into the interaction sphere [1/s]
f_FM=delta*delta*np.pi*c*C	# vapor flux into the interaction sphere in free molecular limit [1/s]
print ("f_FM="+str(f_FM*1e-9)+"[1/ns]\tf_sim="+str(f_sim*1e-9)+"[1/ns]")

os.system("make")

plot=plot.plot()
plot.pltNormal()
plot.fig, plot.axs = plt.subplots(3,2,figsize=(8,10))
for i in np.arange(6):
	plot.axNormal(plot.axs.flat[i])

				## Unknown parameters
#-----------------------------------------------------------------------------#
teq=0.5e-6     # equilibliumed time [s]
tEND=1.5e-6
tcut=1e-13       # cut residence time [s]
Nmax=40         # Max number of sticking vapors
dt_post=1e-11   # dt in analysis, dt_post [s]
delta=1e-8	# diameter of interaction sphere [m]
delta2=(delta*1e10)**2      # square of delta [ang^2]
directory="../../../nucleation/NaCl/Na2Cl_toluene/"
I=10

startTime=2e-9
endTime=5e-9

checkMode=0 # display
figOutput=1

Mvapor=74       # vapor molar mass [g/mol]
M=1/(1/Mvapor)		# vapor mass [kg/mol]
kb=1.38e-23	# boltzmann constant [J/K]
R=8.314		# gas constant	[Jmol/K]
T=300.0		# temperature	[K]
c=(8*R*T/M/np.pi)**0.5	# vapor mean thermal speed [m/s]
pv0=10
#-----------------------------------------------------------------------------#


				## Reading energy file(s)
#-----------------------------------------------------------------------------#
U=np.loadtxt(str(directory)+"U_"+str(I)+".dat")
t_tot=np.max(U.T[0])*1e-15  # total simulation time [s]
negs=np.where((U[:,0]>tEND*1e15))
U=np.delete(U,negs,axis=0)
plot.plotEnergies(U,teq,tEND)
#-----------------------------------------------------------------------------#

				## Reading vapor in out time files
#-----------------------------------------------------------------------------#
inVapor=np.loadtxt(str(directory)+"vapor_in_"+str(I)+".dat")
negs=np.where((inVapor[:,1]>tEND*1e15))
inVapor=np.delete(inVapor,negs,axis=0)

outVapor=np.loadtxt(str(directory)+"vapor_out_"+str(I)+".dat")
negs=np.where((outVapor[:,1]>tEND*1e15))
outVapor=np.delete(outVapor,negs,axis=0)

Npost=int(t_tot/dt_post)					# total number of steps in analysis, Npost=/dt_post

times=np.arange(Npost)*dt_post
Nstick=np.zeros(Npost)	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]
usedLines=np.arange(np.size(outVapor.T[0]))		# Flag for outVapor (equal to -1 if it is used)
ts=np.zeros(np.size(inVapor.T[0]))	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]
tth=np.zeros(np.size(inVapor.T[0]))	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]

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

plot.plotNvap(times,teq,Nstick)
plot.plotStickTimeDist(ts,tth)

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

plot.plotStickVaporDist(nv,ppoi,psim,Nbase)


pv=pv0-Nbase 	# vapor pressure [Pa]
C=pv/kb/T	# vapor concentraiton [1/m3]
f_sim=np.size(ts)/t_tot	# vapor flux into the interaction sphere [1/s]
f_FM=delta*delta*np.pi*c*C	# vapor flux into the interaction sphere in free molecular limit [1/s]
print ("f_FM="+str(f_FM*1e-9)+"[1/ns]\tf_sim="+str(f_sim*1e-9)+"[1/ns]")

#-----------------------------------------------------------------------------#
ionData=np.loadtxt(str(directory)+"ion_300_"+str(I)+".dat")
dtInput=ionData[1][0]-ionData[0][0]
os.system("./mobility.out "+str(I)+" 1 "+str(teq)+" "+str(startTime)+" "+str(endTime)+" "+str(dtInput*1e-9)+" "+str(directory)+" "+str(tEND))
MSDVAF=np.loadtxt(str(directory)+"TIME_MSD_VAF."+str(I))
diffusionData=np.loadtxt(str(directory)+"DiffusionCoefficients."+str(I),skiprows=1)
plot.plotMSDVAF(MSDVAF,diffusionData)
#-----------------------------------------------------------------------------#

plot.plotShow(figOutput,str(directory)+"fig"+str(I)+".png")
