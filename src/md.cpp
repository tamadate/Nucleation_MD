//------------------------------------------------------------------------
#include "md.hpp"
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
	output_initial();
	setCondition(condfile, atomFile);

	d_size=pow(Nof_around_gas*kb*T/p,1/3.0)*1e10;
	V=d_size*d_size*d_size;
	dt = 0.5;	/*	fs	*/
	CUTOFF = 20.0;	/*	A	*/
	MARGIN = 10.0;	/*	A	*/
	ML2 = (CUTOFF+MARGIN)*(CUTOFF+MARGIN);
	CL2 = (CUTOFF*CUTOFF);
	OBSERVE=10000000;


	if(pp->vaportype==1) vars->atomVapor = vars->makeAtomMeOH();
	if(pp->vaportype==2) vars->atomVapor = vars->makeAtomTIP3P();
	vars->read_initial(atomFile);
	vars->set_initial_velocity(pp);
	setPotential(flags);
	initialization_gas();	//Set initial positions & velocities for gas
  initialization_vapor();	//Set initial positions & velocities for vapor
	analysis_ion();
	make_pair();
	margin_length = MARGIN;
	vars->tzero();

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
