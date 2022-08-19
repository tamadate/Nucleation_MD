import numpy as np
import matplotlib.pylab as plt

inVapor=np.loadtxt("vapor_in_1.dat")
outVapor=np.loadtxt("vapor_out_1.dat")

t_tot=0.04e-6	# total simulation time, t_tot [s]
dt=1e-15	# simulation time step, dt [s]
dt_post=1e-11								# dt in analysis, dt_post [s]
Npost=int(t_tot/dt_post)					# total number of steps in analysis, Npost=/dt_post

times=np.arange(0,t_tot,dp_post)
Nstick=np.zeros(Npost)	# Number of vapors [N(t0),N(t0+dt_post),N(t1+2*dt_post),...,N(t_tot)]
usedLines=np.arange(np.size(outVapor.T[0]))		# Flag for outVapor (equal to -1 if it is used)

for i in inVapor:
	id1=i[0]
	time1=int(i[1]/dt_post)
	time2=Npost-1
	for loop in np.arange(np.size(usedLines)):
		if(usedLines[loop]==-1):
			continue
		if(outVapor[loop][0]==id1):
			time2=int(outVapor[loop][1]/dt_post)
			usedLines[loop]=-1
			break
	Nstick[time1:time2]+=1

np.savetxt("stick.dat",Nstick)

plt.plot(dt_post,Nstick)
