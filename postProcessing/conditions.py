import numpy as np

class conditions:
	#-----------------------------------------------------------------------------#
	kb=1.38e-23	# boltzmann constant [J/K]
	R=8.314		# gas constant	[Jmol/K]
	T=300.0		# temperature	[K]

	teq=0.5e-6     # equilibliumed time [s]
	tEND=1.5e-6
	tcut=1e-13       # cut residence time [s]
	Nmax=40         # Max number of sticking vapors
	dt_post=1e-11   # dt in analysis, dt_post [s]
	delta=1e-8	# diameter of interaction sphere [m]
	delta2=(delta*1e10)**2      # square of delta [ang^2]
	directory="../../../nucleation/NaCl/Na2Cl_toluene/"
	#directory="/home/tama3rdgen/vaporUptake/angioII+1/"
	I=9

	startTime=2e-9
	endTime=5e-9

	figOutput=1

	pv0=10
	Mvapor=74       # vapor molar mass [g/mol]
	M=1/(1/Mvapor)		# vapor mass [kg/mol]
	c=(8*R*T/M/np.pi)**0.5	# vapor mean thermal speed [m/s]

	def setMasses(self,mv,mi):
		self.Mvapor=mv
		self.M=1/(1/mv+1/mi)		# vapor mass [kg/mol]
		self.c=(8*R*T/M/np.pi)**0.5	# vapor mean thermal speed [m/s]
	#-----------------------------------------------------------------------------#
