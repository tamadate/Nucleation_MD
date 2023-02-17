import numpy as np
import math
import matplotlib.pyplot as plt
import os
import glob
import plot
from scipy.optimize import minimize

os.system("make")

#-----------------------------------------------------------------------------#
plot=plot.plot()

plot.pltNormal()
plot.fig, plot.axs = plt.subplots(1,1,figsize=(5,5))
plot.axNormal(plot.axs)

plot.labelSize=15
plot.titleSize=18
plot.lineWidth=0.8
#-----------------------------------------------------------------------------#


				## Unknown parameters
#-----------------------------------------------------------------------------#
directories=["/media/tama3rdgen/6TB/vaporUptake/bradykinin/",\
"/media/tama3rdgen/6TB/vaporUptake/angio2+_new/",\
"/media/tama3rdgen/6TB/vaporUptake/angio1+_new/"]

plot.plotAveStickTime(directories)
