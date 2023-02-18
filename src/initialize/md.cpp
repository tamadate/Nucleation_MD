
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
	con = new MDcondition();
	obs = new Observer(vars,con,calcNumber);
	pp = new Physical();
	flags = new FLAG();
	thermo = new Thermostat();

	dt = 0.5;	/*	fs	*/

	readCondFile(condfile);
	pp->readIonProp(atomFile);
	pp->readVaporProp(vaporFile);
	pp->setPhysicalProp(gastype);

	vars->readIonFile(atomFile);
	vars->readVaporFile(vaporFile);
	vars->setBMHPotential();
	vars->setCrossPotentials();
	vars->ionRotation();
	vars->ionInitialVelocity(pp->T);
	for (auto &a : InterInter) a->initial(vars);
	for (auto &a : IntraInter) a->initial(vars);
	for (auto &a : IntraInter) a->printName();
	for (auto &a : InterInter) a->printName();

	initialization_gas();	//Set initial positions & velocities for gas
  	initializatIonVapor();	//Set initial positions & velocities for vapor
	getIonCenterProp();
	make_pair();
	vars->tzero();

	mbdist = new MBdist(pp->cgas,pp->mgas,pp->T);
	mbdistV = new MBdist(pp->cvapor,pp->mvapor,pp->T);
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

