import numpy as np
import os
import plot
import diffusionCoeff
import conditions
import cluster
import stickPosition
from scipy.optimize import minimize

con=conditions.conditions

con.teq=0.5e-6     # equilibliumed time [s]
con.tEND=1.5e-6
con.tcut=1e-13       # cut residence time [s]
con.Nmax=40         # Max number of sticking vapors
con.dt_post=1e-11   # dt in analysis, dt_post [s]
con.startTime=2e-9
con.endTime=5e-9

con.directory="/home/tama3rdgen/vaporUptake/angioII+2/"
#con.directory="../../../nucleation/valine/"
con.I=10
con.figOutput=0

#clu=cluster.cluster(con)
#clu.compute()
#diff=diffusionCoeff.diffusionCoeff(con)
#diff.compute()
stk=stickPosition.stickPosition(con)
stk.compute()
