import numpy as np
import os
import plot
import diffusionCoeff
import conditions
import cluster
import stickPosition
from scipy.optimize import minimize

con=conditions.conditions

con.teq=0.0e-6     # equilibliumed time [s]
con.tEND=0.15e-6
con.tcut=1e-9       # cut residence time [s]
con.Nmax=20         # Max number of sticking vapors
con.dt_post=1e-11   # dt in analysis, dt_post [s]
con.startTime=5e-9
con.endTime=10e-9

#con.directory="/home/tama3rdgen/vaporUptake/angioII+2/"
con.directory="/home/tama3rdgen/test2/valine/500Pa_MeOH/"
con.directory="/home/tama3rdgen/vaporUptake/bradykinin+2/1/"

con.I=2
con.figOutput=0
con.pv=500
con.setMasses(con,32e-3,117.15e-3)
con.calcParams(con)


clu=cluster.cluster(con)
clu.Upot()
clu.vapor()
'''
clu=cluster.cluster(con)
clu.compute()
'''
diff=diffusionCoeff.diffusionCoeff(con)
diff.compute()
#stk=stickPosition.stickPosition(con)
#stk.compute()
