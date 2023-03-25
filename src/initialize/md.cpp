
#include "../md.hpp"


/////////////////////////////////////////////////////////////////////
/*
	constructor
*/
/////////////////////////////////////////////////////////////////////
MD::MD(char* condfile, int calcNumber) {
	vars = new Variables();
	vars->calcID =  calcNumber;
	con = new MDcondition();
	obs = new Observer(vars,con);
	pp = new Physical();

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
	for (auto &a : funcs) a->initial(vars, pp, con, obs);
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
	delete mbdist;
	delete mbdistV;
	for (auto &a : funcs) delete a;
}

