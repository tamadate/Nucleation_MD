//------------------------------------------------------------------------
#include "../md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*
	constructor
*/
/////////////////////////////////////////////////////////////////////
MD::MD(char* condfile, int calcNumber) {
	startTime=omp_get_wtime();
	calculation_number =  calcNumber;
	#pragma omp parallel
	{
		#pragma omp single
		{
			Nth=omp_get_num_threads();
		}
	}
	vars = new Variables();
	obs = new Observer();
	pp = new Physical();
	flags = new FLAG();
	mbdist = new MBdist();
	mbdistV = new MBdist();

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
	output_initial();

	vars->readIonFile(atomFile);
	vars->readVaporFile(vaporFile);
	vars->setBMHPotential();
	vars->setCrossPotentials();
	vars->ionRotation();
	vars->ionInitialVelocity(T);

	setPotential(flags,1);
	initialization_gas();	//Set initial positions & velocities for gas
  	initialization_vapor();	//Set initial positions & velocities for vapor
	analysis_ion();
	make_pair();
	margin_length = MARGIN;
	vars->tzero();

	mbdist -> makeWeightedMB(pp->cgas,pp->mgas,T);
	mbdistV -> makeWeightedMB(pp->cvapor,pp->mvapor,T);
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


/*########################################################################################

-----Calculating center of mass and velocity of center of mass-----

#######################################################################################*/

void
MD::analysis_ion(void) {
	for ( int i = 0; i < 3; i++ ) ion_r[ i ] = 0, ion_v[ i ] = 0;
	for ( auto &a : vars -> ions)	{
		ion_r[ 0 ] += a.qx * a.mass;
		ion_r[ 1 ] += a.qy * a.mass;
		ion_r[ 2 ] += a.qz * a.mass;
		ion_v[ 0 ] += a.px * a.mass;
		ion_v[ 1 ] += a.py * a.mass;
		ion_v[ 2 ] += a.pz * a.mass;
	}
	for ( int i = 0; i < 3; i++ ) ion_r[ i ] /= pp -> Mion;
	for ( int i = 0; i < 3; i++ ) ion_v[ i ] /= pp -> Mion;
}

void
MD::analysis_gas(void) {
	for ( int i = 0; i < 3; i++) gas_r[ i ] = 0, gas_v[ i ] = 0;
	for (auto &a : vars -> gases)	{
		gas_r[ 0 ] += a.qx * a.mass;
		gas_r[ 1 ] += a.qy * a.mass;
		gas_r[ 2 ] += a.qz * a.mass;
		gas_v[ 0 ] += a.px *a.mass;
		gas_v[ 1 ] += a.py *a.mass;
		gas_v[ 2 ] += a.pz *a.mass;
	}
	for ( int i = 0; i < 3; i++ ) gas_r[ i ] /= pp -> Mgas;
	for ( int i = 0; i < 3; i++) gas_v[ i ] /= pp -> Mgas;
}
