import numpy as np
import math
import matplotlib.pyplot as plt
import os
from scipy.optimize import minimize


class plot:

    labelSize=10
    titleSize=12
    lineWidth=0.8

    #-----------------------------------------------------------------------------#
    def pltNormal(self):
        plt.rcParams['ytick.direction'] = 'in'
        plt.rcParams['xtick.direction'] = 'in'
        #plt.rcParams['figure.subplot.bottom'] = 0.2
        #plt.rcParams['figure.subplot.left'] = 0.2
        #plt.rcParams['font.family'] = 'Arial'
        plt.rcParams["font.size"]=8

    def axNormal(self,ax):
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.tick_params(axis='x')
        ax.tick_params(axis='y')

    				## Reading energy file(s)
    #-----------------------------------------------------------------------------#

    def plotEnergies(self,U,teq,tEND):
        #K=np.loadtxt("K_1.dat")
        labels=["Uion","Ugas","Uvap","Ugi","Ugg","Uvg","Uvi","Uvv"]
        display=[1,0,0,0,0,0,1,1]
        colors=["blue","blue","blue","blue","blue","blue","aqua","navy"]

        t_tot=tEND  # total simulation time [s]
        Xmax=int(t_tot/200e-9)*200
        Xmax*=1.4

        #axs.flat[0].set_xlim([0,Xmax])
        self.axs.flat[0].set_title("(a) Potential Energies",loc='left',fontsize=self.titleSize)
        self.axs.flat[0].set_xlabel("Time [ns]",fontsize=self.labelSize)
        self.axs.flat[0].set_ylabel("Energy [kcal/mol]",fontsize=self.labelSize)
        for i in np.arange(np.size(U[0])-1):
            if(display[i]):
        	       self.axs.flat[0].plot(U.T[0]*1e-6,U.T[i+1],label=labels[i],linewidth=self.lineWidth,color=colors[i])
        self.axs.flat[0].plot(U.T[0]*1e-6,np.sum(U.T[1:],axis=0),label="Total",linewidth=self.lineWidth,color="black")
        self.axs.flat[0].axvline(x = teq*1e9, ls='--', color = 'black',linewidth=self.lineWidth)
        self.axs.flat[0].legend()

    				## Reading vapor in out time files
    #-----------------------------------------------------------------------------#
    def plotNvap(self,times,teq,Nstick):
        self.axs.flat[1].set_title("(b) MeOHs in effective domain",loc='left',fontsize=self.titleSize)
        self.axs.flat[1].set_xlabel("Time [ns]",fontsize=self.labelSize)
        self.axs.flat[1].set_ylabel("$\it {N}$$_ {vap}$ [-]",fontsize=self.labelSize)
        self.axs.flat[1].axvline(x = teq*1e9, color = 'black', ls="--",linewidth=self.lineWidth)
        self.axs.flat[1].scatter(times*1e9,Nstick,color="black",s=20)

    #-----------------------------------------------------------------------------#
    def plotStickTimeDist(self,ts,tth):
        negs=np.where(ts<=0)
        self.axs.flat[2].set_title("(c) Vapor sticking time distribution",loc='left',fontsize=self.titleSize)
        self.axs.flat[2].set_xlabel("Logarithm of time [-]",fontsize=self.labelSize)
        self.axs.flat[2].set_ylabel("Number of event [-]",fontsize=self.labelSize)
        self.axs.flat[2].set_yscale("log")
        self.axs.flat[2].hist(np.log10(np.delete(ts,negs)),alpha=0.3,bins=50,label="$\it {t}$$_{sim}$",color="blue")
        self.axs.flat[2].hist(np.log10(np.delete(tth,negs)),alpha=0.3,bins=30,label="$\it {t}$$_ {sim}$-$\it t$$_ {th}$",color="black")
        #axs.flat[2].axvline(x = np.log10(tcut), color = 'black', ls="--",linewidth=lineWidth)
        self.axs.flat[2].legend()

    #-----------------------------------------------------------------------------#
    def plotStickVaporDist(self,nv,ppoi,psim,Nbase):
        self.axs.flat[3].set_title("(d) Sticking vapor distribution",loc='left',fontsize=self.titleSize)
        self.axs.flat[3].set_title("")
        self.axs.flat[3].set_xlabel("Number of sticking vapor",fontsize=self.labelSize)
        self.axs.flat[3].set_ylabel("Frequency [-]",fontsize=self.labelSize)
        self.axs.flat[3].set_yscale("log")
        self.axs.flat[3].set_ylim([1e-5,1])
        self.axs.flat[3].scatter(nv,ppoi,label="Poisson",s=20,color="black")
        self.axs.flat[3].scatter(nv,psim,label="Simulation",s=20,marker="v",color="blue")
        self.axs.flat[3].scatter(nv-Nbase,psim,label="Simulation("+str(int(-Nbase))+")",s=20,marker="^",color="cyan")
        self.axs.flat[3].legend()

    #-----------------------------------------------------------------------------#
    def plotMSDVAF(self,MSDVAF,diffusionData):
        self.axs.flat[4].set_title("(e) Normalized MSD",loc='left',fontsize=self.titleSize)
        self.axs.flat[4].set_xlabel("Time [ns]",fontsize=self.labelSize)
        self.axs.flat[4].set_ylabel("Normalized MSD [-]",fontsize=self.labelSize)
        self.axs.flat[4].scatter(MSDVAF.T[0],MSDVAF.T[1]/np.max(MSDVAF.T[1]),color="black",s=10)
        self.axs.flat[4].text(1, 0.8, "K = "+'{:.2f}'.format(diffusionData[3]*diffusionData[5])+" cm$^ 2$/Vs", fontsize = 10)
        self.axs.flat[5].set_title("(f) Normalized VAF",loc='left',fontsize=self.titleSize)
        self.axs.flat[5].set_xlabel("Time [ns]",fontsize=self.labelSize)
        self.axs.flat[5].set_ylabel("Normalized VAF [-]",fontsize=self.labelSize)
        self.axs.flat[5].scatter(MSDVAF.T[0],MSDVAF.T[2]/MSDVAF.T[2][0],color="black",s=10)
        self.axs.flat[5].text(2, 0.8, "K = "+'{:.2f}'.format(diffusionData[4]*diffusionData[5])+" cm$^ 2$/Vs", fontsize = 10)

    #-----------------------------------------------------------------------------#
    def plotShow(self,figOutput,fileName):
        self.fig.tight_layout()
        if(figOutput):
            plt.savefig(fileName, dpi=1000)
        plt.show()


    def plotMobilityShift(self,datas,directory):
        #axs.flat[0].set_xlim([0,Xmax])
        self.axs.flat[0].set_title("(a) Inverse of mobility shift",loc='left',fontsize=self.titleSize)
        self.axs.flat[0].set_xlabel("Vapor pressure [Pa]",fontsize=self.labelSize)
        self.axs.flat[0].set_ylabel(r"Mobility [cm$^2$/Vs]",fontsize=self.labelSize)
        for data in datas:
            self.axs.flat[0].scatter(data[0]*10,(data[4]*data[6])**-1,color = "black")

        self.axs.flat[1].set_title("(b) Mobility shift",loc='left',fontsize=self.titleSize)
        self.axs.flat[1].set_xlabel("Vapor pressure [Pa]",fontsize=self.labelSize)
        self.axs.flat[1].set_ylabel(r"Normalized mobility, $Z_p/Z_{p,0}$ [-]",fontsize=self.labelSize)
        for data in datas:
            self.axs.flat[1].scatter(data[0]*10,(data[4]/datas[0][4]),color = "black")

        plt.savefig(str(directory)+"diffusionSummary.png", dpi=1000)
        plt.show()
