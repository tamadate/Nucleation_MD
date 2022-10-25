import numpy as np
import math
import matplotlib.pyplot as plt
import os
from scipy.optimize import minimize


class plot:

    labelSize=12
    titleSize=15
    lineWidth=1
    fs=(5,5)

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

    def plotEnergies(self,U,teq,tEND,figOutput):
        self.pltNormal()
        fig, axs = plt.subplots(1,1,figsize=self.fs)
        self.axNormal(axs)

        #K=np.loadtxt("K_1.dat")
        labels=["Uion","Ugas","Uvap","Ugi","Ugg","Uvg","Uvi","Uvv"]
        display=[1,0,0,0,0,0,1,1]
        colors=["blue","blue","blue","blue","blue","blue","aqua","navy"]

        t_tot=tEND  # total simulation time [s]
        Xmax=int(t_tot/200e-9)*200
        Xmax*=1.4

        #axs.set_xlim([0,Xmax])
        axs.set_xlabel("Time [ns]",fontsize=self.labelSize)
        axs.set_ylabel("Energy [kcal/mol]",fontsize=self.labelSize)
        for i in np.arange(np.size(U[0])-1):
            if(display[i]):
        	       axs.plot(U.T[0]*1e-6,U.T[i+1],label=labels[i],linewidth=self.lineWidth,color=colors[i])
        axs.plot(U.T[0]*1e-6,np.sum(U.T[1:],axis=0),label="Total",linewidth=self.lineWidth,color="black")
        axs.axvline(x = teq*1e9, ls='--', color = 'black',linewidth=self.lineWidth)
        axs.legend()
        fig.tight_layout()
        if(figOutput):
            plt.savefig(fileName, dpi=1000)
        plt.show()

    				## Reading vapor in out time files
    #-----------------------------------------------------------------------------#
    def plotNvap(self,times,teq,Nstick,figOutput):
        self.pltNormal()
        fig, axs = plt.subplots(1,1,figsize=(5,5))
        self.axNormal(axs)
        axs.set_xlabel("Time [ns]",fontsize=self.labelSize)
        axs.set_ylabel("$\it {N}$$_ {vap}$ [-]",fontsize=self.labelSize)
        axs.axvline(x = teq*1e9, color = 'black', ls="--",linewidth=self.lineWidth)
        axs.scatter(times*1e9,Nstick,color="black",s=20)
        fig.tight_layout()
        if(figOutput):
            plt.savefig(fileName, dpi=1000)
        plt.show()

    #-----------------------------------------------------------------------------#
    def plotStickTimeDist(self,ts,tth,figOutput):
        self.pltNormal()
        fig, axs = plt.subplots(1,1,figsize=(5,5))
        self.axNormal(axs)
        negs=np.where(ts<=0)
        axs.set_xlabel("Logarithm of time [-]",fontsize=self.labelSize)
        axs.set_ylabel("Number of event [-]",fontsize=self.labelSize)
        axs.set_yscale("log")
        axs.hist(np.log10(np.delete(ts,negs)),alpha=0.3,bins=50,label="$\it {t}$$_{sim}$",color="blue")
        axs.hist(np.log10(np.delete(tth,negs)),alpha=0.3,bins=30,label="$\it {t}$$_ {sim}$-$\it t$$_ {th}$",color="black")
        #axs.axvline(x = np.log10(tcut), color = 'black', ls="--",linewidth=lineWidth)
        axs.legend()
        fig.tight_layout()
        if(figOutput):
            plt.savefig(fileName, dpi=1000)
        plt.show()

    #-----------------------------------------------------------------------------#
    def plotStickVaporDist(self,nv,ppoi,psim,Nbase,figOutput):
        self.pltNormal()
        fig, axs = plt.subplots(1,1,figsize=(5,5))
        self.axNormal(axs)
        axs.set_xlabel("Number of sticking vapor",fontsize=self.labelSize)
        axs.set_ylabel("Frequency [-]",fontsize=self.labelSize)
        axs.set_yscale("log")
        axs.set_ylim([1e-5,1])
        axs.scatter(nv,ppoi,label="Poisson",color="black")
        axs.scatter(nv,psim,label="Simulation",color="blue")
        #axs.scatter(nv-Nbase,psim,label="Simulation("+str(int(-Nbase))+")",marker="^",color="cyan")
        axs.legend()
        fig.tight_layout()
        if(figOutput):
            plt.savefig(fileName, dpi=1000)
        plt.show()

    #-----------------------------------------------------------------------------#
    def plotMSDVAF(self,MSDVAF,diffusionData,figOutput):
        self.pltNormal()
        fig, axs = plt.subplots(1,1,figsize=(5,5))
        self.axNormal(axs)
        axs.set_title("(e) Normalized MSD",loc='left',fontsize=self.titleSize)
        axs.set_xlabel("Time [ns]",fontsize=self.labelSize)
        axs.set_ylabel("Normalized MSD [-]",fontsize=self.labelSize)
        axs.scatter(MSDVAF.T[0],MSDVAF.T[1]/np.max(MSDVAF.T[1]),color="black",s=10)
        axs.text(1, 0.8, "K = "+'{:.2f}'.format(diffusionData[3]*diffusionData[5])+" cm$^ 2$/Vs", fontsize = 10)
        axs.set_title("(f) Normalized VAF",loc='left',fontsize=self.titleSize)
        axs.set_xlabel("Time [ns]",fontsize=self.labelSize)
        axs.set_ylabel("Normalized VAF [-]",fontsize=self.labelSize)
        axs.scatter(MSDVAF.T[0],MSDVAF.T[2]/MSDVAF.T[2][0],color="black",s=10)
        axs.text(2, 0.8, "K = "+'{:.2f}'.format(diffusionData[4]*diffusionData[5])+" cm$^ 2$/Vs", fontsize = 10)
        fig.tight_layout()
        if(figOutput):
            plt.savefig(fileName, dpi=1000)
        plt.show()

    #-----------------------------------------------------------------------------#
    def plotShow(self,figOutput,fileName):
        fig.tight_layout()
        if(figOutput):
            plt.savefig(fileName, dpi=1000)
        plt.show()

    def plotMobilityShift(self,datas,directory):
        self.pltNormal()
        fig, axs = plt.subplots(1,1,figsize=(5,5))
        self.axNormal(axs)
        #axs.set_xlim([0,Xmax])
        axs.set_title("(a) Inverse of mobility shift",loc='left',fontsize=self.titleSize)
        axs.set_xlabel("Vapor pressure [Pa]",fontsize=self.labelSize)
        axs.set_ylabel(r"Mobility [cm$^2$/Vs]",fontsize=self.labelSize)
        for data in datas:
            self.scatter(data[0]*10,(data[4]*data[6])**-1,color = "black")

        axs.set_title("(b) Mobility shift",loc='left',fontsize=self.titleSize)
        axs.set_xlabel("Vapor pressure [Pa]",fontsize=self.labelSize)
        axs.set_ylabel(r"Normalized mobility, $Z_p/Z_{p,0}$ [-]",fontsize=self.labelSize)
        for data in datas:
            axs.scatter(data[0]*10,(data[4]/datas[0][4]),color = "black")

        plt.savefig(str(directory)+"diffusionSummary.png", dpi=1000)
        plt.show()
