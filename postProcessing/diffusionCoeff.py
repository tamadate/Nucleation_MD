import numpy as np
import matplotlib.pyplot as plt
import os
import plot
import conditions



class diffusionCoeff:
    def __init__(self,con):
        self.plot=plot.plot()
        self.plot.pltNormal()
        self.plot.directory=con.directory
        self.con=con

    def compute(self):
        os.system("make")
        ionData=np.loadtxt(self.con.directory+"ion_300_"+str(self.con.I)+".dat")
        dtInput=ionData[1][0]-ionData[0][0]
        os.system("./mobility.out "+str(self.con.I)+" 1 "+str(self.con.teq)+" "+str(self.con.startTime)+" "+\
        str(self.con.endTime)+" "+str(dtInput*1e-9)+" "+self.con.directory+" "+str(self.con.tEND))
        MSDVAF=np.loadtxt(self.con.directory+"TIME_MSD_VAF."+str(self.con.I))
        diffusionData=np.loadtxt(self.con.directory+"DiffusionCoefficients."+str(self.con.I),skiprows=1)

        self.plot.plotMSDVAF(MSDVAF,diffusionData,self.con.figOutput)
        #self.plot.plotShow(self.con.figOutput,self.con.directory+"fig"+str(self.con.I)+".png")
