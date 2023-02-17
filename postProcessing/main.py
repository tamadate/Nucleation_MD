import numpy as np
import os
import plot
import diffusionCoeff
import conditions
import cluster
import error
import stickPosition
from scipy.optimize import minimize

con=conditions.conditions

con.teq=0.0e-6     # equilibliumed time [s]
con.tEND=0.23e-6
con.tcut=1e-13       # cut residence time [s]
con.Nmax=11         # Max number of sticking vapors
con.dt_post=1e-11   # dt in analysis, dt_post [s]
con.startTime=5e-9
con.endTime=10e-9

#con.directory="/home/tama3rdgen/vaporUptake/angioII+2/"
con.directory="/media/tama3rdgen/6TB/vaporUptake/bradykinin/v10/100Pa/"
#con.directory="/media/tama3rdgen/6TB/vaporUptake/angio2+_new/v14/140Pa/"
#con.directory="/media/tama3rdgen/6TB/vaporUptake/angio1+_new/v8/80Pa/"
print(con.directory)
con.I=1
con.figOutput=1
con.pv0=100
con.pal=100
con.setMasses(con,32e-3,117.15e-3)
con.calcParams(con)

error.errorCheck(con)

clu=cluster.cluster(con)
#clu.compute()
clu.Upot()
clu.evaporationRate()
#print(clu.con.error)

#diff=diffusionCoeff.diffusionCoeff(con)
#diff.compute()

#stk=stickPosition.stickPosition(con)
#stk.compute()
