
#include "../md.hpp"

MD::MD(char* condfile, int calcNumber) {
	vars = new Variables();		// generate class for variable 
	vars->calcID =  calcNumber;	// setting calculation ID
	con = new MDcondition();	// generate class calculation conditions
	obs = new Observer(vars,con);	// generate class for observer
	pp = new Physical();			// generate class for physical properties 

	readCondFile(condfile);	// read conditions file

	pp->readIonProp(atomFile);	// read ion file and set physical properties
	pp->readVaporProp(vaporFile);	// read vapor file and set physical properties
	pp->setPhysicalProp(gastype);	// set all physical properties

	// read ion and vapor files, storing molecular structures and potential parameters
	vars->readIonFile(atomFile);	
	vars->readVaporFile(vaporFile);	
	vars->setCrossPotentials();	// set cross (gas-ion, gas-vapor, and ion-vapor) potential parameters
	vars->ionRotation();		// rotate ion
	vars->ionInitialVelocity(pp->T);	// set initial velocites for all atom in the ion

	// initializing potential functions and print the potential names used in simulation
	for (auto &a : InterInter) a->initial(vars);
	for (auto &a : IntraInter) a->initial(vars);
	for (auto &a : IntraInter) a->printName();
	for (auto &a : InterInter) a->printName();

	initialization_gas();	// Set initial positions & velocities for gas
  	initializatIonVapor();	// Set initial positions & velocities for vapor
	getIonCenterProp();		// get ion COM and linear velocity
	make_pair();			// make pair list with initial positions

	// generate MB distribution class
	mbdist = new MBdist(pp->mgas,pp->T);	
	mbdistV = new MBdist(pp->mvapor,pp->T);

	// initializing all functions
	for (auto &a : funcs) a->initial(vars, pp, con, obs);
}


MD::~MD(void) {
	delete vars;
	delete obs;
	delete pp;
	delete mbdist;
	delete mbdistV;
	for (auto &a : funcs) delete a;
}

