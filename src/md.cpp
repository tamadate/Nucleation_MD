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

	dt = 0.5;	/*	fs	*/
	CUTOFF = 20.0;	/*	A	*/
	MARGIN = 10.0;	/*	A	*/
	OBSERVE=10000000;
	T=300;
	p=1e5;
	positionLogStep=0;
	totalVaporIn=0;
	setCondition(condfile);
	output_initial();

	vars->read_initial(atomFile,vaporFile);
	vars->ionInitialVelocity(T);

	setPotential(flags,1);
	if(flags->takeOver==1){
		takeOver();
		//cout<<takeOverFile<<endl;
	}
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
	Molecule *ion = vars -> effectiveIn[0].data();
	ion[0].qx=0;
	ion[0].qy=0;
	ion[0].qz=0;
	ion[0].px=0;
	ion[0].py=0;
	ion[0].pz=0;

	for ( auto &a : ion[0].inAtoms)	{
		ion[0].qx += a.qx * a.mass;
		ion[0].qy += a.qy * a.mass;
		ion[0].qz += a.qz * a.mass;
		ion[0].px += a.px * a.mass;
		ion[0].py += a.py * a.mass;
		ion[0].pz += a.pz * a.mass;
	}
	ion[0].qx /= pp -> Mion;
	ion[0].qy /= pp -> Mion;
	ion[0].qz /= pp -> Mion;
	ion[0].px /= pp -> Mion;
	ion[0].py /= pp -> Mion;
	ion[0].pz /= pp -> Mion;
}
