import numpy as np
import math
import os
import plot
import conditions
import matplotlib.pyplot as plt
from scipy.optimize import minimize


def errorCheck(con):
	con.error=[]
	f=open(con.directory+"error.dat","w")
	for i in np.arange(con.pal):
		U=np.loadtxt(con.directory+"U_"+str(con.I+i)+".dat")
		if(np.count_nonzero(np.isnan(U))>0):
			con.error=np.append(con.error,int(i))
			#con.pal-=1
			f.write(str(int(i))+"\n")
		if(np.max(U.T[0])<con.tEND*1e15):
			con.error=np.append(con.error,int(i))
			#con.pal-=1
			f.write(str(int(i))+"\n")


def errorAdd(con,neg,indexes):
	f=open(con.directory+"error.dat","a")
	for n in neg:
		con.error=np.append(con.error,indexes[n])
		#con.pal-=1
		f.write(str(int(indexes[n]))+"\n")
