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
        plt.rcParams['figure.subplot.bottom'] = 0.2
        plt.rcParams['figure.subplot.left'] = 0.2
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
            plt.savefig(str(self.directory)+"Energy.png", dpi=1000)
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
            plt.savefig(str(self.directory)+"Nvap.png", dpi=1000)
        plt.show()

    #-----------------------------------------------------------------------------#
    def plotStickTimeDist(self,tss,tth,figOutput):
        self.pltNormal()
        fig, axs = plt.subplots(1,1,figsize=(5,5))
        self.axNormal(axs)
        negs=np.where(tss<=0)
        axs.set_xlabel("Logarithm of time [-]",fontsize=self.labelSize)
        axs.set_ylabel("Number of event [-]",fontsize=self.labelSize)
        #axs.set_yscale("log")
        axs.hist(np.log10(np.delete(tss,negs)),alpha=0.3,bins=50,label="$\it {t}$$_{sim}$",color="blue")
        ymin,ymax=axs.get_ylim()
        axs.text(-10.5, ymax*0.8, r"$t_s$ = "+'{:.3f}'.format(np.average(np.delete(tss,negs))*1e9)+" ns", fontsize = 10)
        #axs.hist(np.log10(np.delete(tth,negs)),alpha=0.3,bins=30,label="$\it {t}$$_ {sim}$-$\it t$$_ {th}$",color="black")
        #axs.axvline(x = np.log10(tcut), color = 'black', ls="--",linewidth=lineWidth)
        #axs.legend()
        fig.tight_layout()
        if(figOutput):
            plt.savefig(str(self.directory)+"stickTimeDist.png", dpi=1000)
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
        print("Poisson")
        print(ppoi)
        print("Simulation")
        print(psim)
        #axs.scatter(nv-Nbase,psim,label="Simulation("+str(int(-Nbase))+")",marker="^",color="cyan")
        axs.legend()
        fig.tight_layout()
        if(figOutput):
            plt.savefig(str(self.directory)+"vaporDist.png", dpi=1000)
        plt.show()

    #-----------------------------------------------------------------------------#
    def plotMSDVAF(self,MSDVAF,diffusionData,figOutput):
        self.pltNormal()
        fig, axs = plt.subplots(1,2,figsize=(10,5))
        for ax in axs.flat:
            self.axNormal(ax)
        axs.flat[0].set_xlabel("Time [ns]",fontsize=self.labelSize)
        axs.flat[0].set_ylabel("Normalized MSD [-]",fontsize=self.labelSize)
        axs.flat[0].scatter(MSDVAF.T[0],MSDVAF.T[1]/np.max(MSDVAF.T[1]),color="black",s=10)
        axs.flat[0].text(1, 0.8, "K = "+'{:.2f}'.format(diffusionData[3]*diffusionData[5])+" cm$^ 2$/Vs", fontsize = 10)
        axs.flat[1].set_xlabel("Time [ns]",fontsize=self.labelSize)
        axs.flat[1].set_ylabel("Normalized VAF [-]",fontsize=self.labelSize)
        axs.flat[1].scatter(MSDVAF.T[0],MSDVAF.T[2]/MSDVAF.T[2][0],color="black",s=10)
        axs.flat[1].text(2, 0.8, "K = "+'{:.2f}'.format(diffusionData[4]*diffusionData[5])+" cm$^ 2$/Vs", fontsize = 10)
        fig.tight_layout()
        if(figOutput):
            plt.savefig(str(self.directory)+"MSD_VAF.png", dpi=1000)
        plt.show()

    #-----------------------------------------------------------------------------#
    def plotShow(self,figOutput,fileName):
        fig.tight_layout()
        if(figOutput):
            plt.savefig(fileName, dpi=1000)
        plt.show()

    def plotMobilityShift(self,press,datas,directory,exp):
        self.axs.set_xlabel("Vapor pressure [Pa]",fontsize=self.labelSize)
        #self.axs.set_ylabel(r"Normalized mobility, $Z_p/Z_{p,0}$ [-]",fontsize=self.labelSize)
        self.axs.set_ylabel(r"Mobility shift, $Z_p/Z_{p,0}$ [-]",fontsize=self.labelSize)
        #self.axs.scatter(press,datas.T[3]*datas.T[5],color = "black")
        #self.axs.scatter(press,datas.T[4]*datas.T[5],color = "red")
        self.axs.scatter(press,datas.T[3]*2/(datas[0][3]+datas[0][4]),color = "black")
        self.axs.scatter(press,datas.T[4]*2/(datas[0][3]+datas[0][4]),color = "red")
        self.axs.scatter(exp[0],exp[1],facecolor="none",edgecolor="blue",marker="^")
        self.axs.set_xlim(-10,300)

        plt.savefig(str(directory)+"diffusionSummary.png", dpi=1000)
        plt.show()

    def plotStickLocation(self,vaporIDstr,dist,aminoNames,distAmino,directory):
        self.axs.flat[0].bar(vaporIDstr,dist/np.sum(dist))
        self.axs.flat[0].set_xticklabels(vaporIDstr, rotation = 45)
        self.axs.flat[0].set_xlabel("Atom index",fontsize=self.labelSize)
        self.axs.flat[0].set_ylabel("Probability",fontsize=self.labelSize)
        self.axs.flat[1].bar(aminoNames,distAmino/np.sum(distAmino))
        self.axs.flat[1].set_xlabel("Amino acids name",fontsize=self.labelSize)
        self.axs.flat[1].set_ylabel("Probability",fontsize=self.labelSize)

        plt.savefig(str(directory)+"stickLocation.png", dpi=1000)
        plt.show()
