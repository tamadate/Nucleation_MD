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
teq=0.2e-6     # equilibliumed time [s]
dt=1e-15	    # simulation time step, dt [s]
tsample=1e-9    # time of sample
tshift=1e-12    # initial time shift
directory="../../../test/"
I=1
#-----------------------------------------------------------------------------#



				## Reading energy file(s)
#-----------------------------------------------------------------------------#

i=1
data=np.loadtxt(str(directory)+"ion_300_"+str(i)+".dat")
dt_data=(data[1][0]-data[0][0])*1e-6    # dt of data array [s]
shift_step=int(dt_data/tshift)
t_max=np.max(data.T[0])*1e-9            # maximum time [s]
data=data[int(np.round(teq/dt_data)):]  # delete before equilibrium
Nsample=int(t_max/tshift)-int(tshift/tsample)              # number of samples
print(Nsample)
N_one_sample=int(tsample/dt_data)      # number of data in a sample
MSD=np.zeros(N_one_sample)

for n in np.arange(Nsample):
    print(data[n*shift_step])
    print(n)
