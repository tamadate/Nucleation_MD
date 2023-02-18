
#include "../md.hpp"


/////////////////////////////////////////////////////////////////////
/*
	constructor
*/
/////////////////////////////////////////////////////////////////////
MD::MD(char* condfile, int calcNumber) {
	calculation_number =  calcNumber;
	#pragma omp parallel
	{
		#pragma omp single
		{
			Nth=omp_get_num_threads();
		}
	}
	vars = new Variables();
	obs = new Observer(vars,calcNumber);
	pp = new Physical();
	flags = new FLAG();

	dt = 0.5;	/*	fs	*/
	CUTOFF = 20.0;	/*	A	*/
	MARGIN = 10.0;	/*	A	*/
	OBSERVE=10000000;
	T=300;
	p=1e5;
	positionLogStep=0;


	readCondFile(condfile);
	pp->readIonProp(atomFile);
	pp->readVaporProp(vaporFile);
	pp->setPhysicalProp(gastype,T,p);

	vars->readIonFile(atomFile);
	vars->readVaporFile(vaporFile);
	vars->setBMHPotential();
	vars->setCrossPotentials();
	vars->ionRotation();
	vars->ionInitialVelocity(T);
	for (auto &a : InterInter) a->initial(vars);
	for (auto &a : IntraInter) a->initial(vars);

	initialization_gas();	//Set initial positions & velocities for gas
  	initialization_vapor();	//Set initial positions & velocities for vapor
	getIonCenterProp();
	make_pair();
	margin_length = MARGIN;
	vars->tzero();

	mbdist = new MBdist(pp->cgas,pp->mgas,T);
	mbdistV = new MBdist(pp->cvapor,pp->mvapor,T);
}

/////////////////////////////////////////////////////////////////////
/*
	destructor
*/
/////////////////////////////////////////////////////////////////////
MD::~MD(void) {
	delete vars;
	delete obs;
	delete pp;
	delete flags;
	delete mbdist;
	delete mbdistV;
}

