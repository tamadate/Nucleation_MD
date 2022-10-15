#include "md.hpp"

/*########################################################################################

-----Initialization-----
intialization:
This intialize the calculation. Reading initial positions and setting initial velocities.
reintialization:
This makes connection between thermal relaxation and main diffusion coeficient calculation. It reset position, time, pair list and margine length.

#######################################################################################*/

/////////////////////////////////////////////////////////////////////
/*
	- Randomly arraying gas molecule around an ion with avoiding the
	overlapping. The velocity is picked from  the Maxwell-Boltzumann
	distribution
	- Set ion's center of mass (maybe -> 0), make pair list for initial
	step of simulation, reset the margine size.
*/
/////////////////////////////////////////////////////////////////////

void
MD::initialization_gas(void) {
	vars->effectiveOut[1].resize(Nof_around_gas);
	double nion=vars->effectiveIn[0][0].inAtoms.size();
	double dx,dy,dz, d;
	double dis=5;	/*	minimum gas-gas, gas-ion distance */

  // Set random number generator
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-d_size*0.5,d_size*0.5);

    // Generate random x, y, z positions and calculate minimum gas-gas distance.
	int i=0;
	do {
		Molecule_out a;
		double min_dis=10000.0;
		a.qx=r(mt), a.qy=r(mt), a.qz=r(mt);
		if(i>0){
			for(auto &b : vars->effectiveOut[1]){
				dx=a.qx-b.qx;
				dy=a.qy-b.qx;
				dz=a.qx-b.qx;
				adjust_periodic(dx, dy, dz, d_size);
				d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<min_dis) min_dis=d; // minimum gas-gas distance
			}
		}
		for(auto &b : vars->effectiveIn[0][0].inAtoms) {
			dx=a.qx-b.qx;
			dy=a.qy-b.qy;
			dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz, d_size);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_dis) min_dis=d; // minimum gas-ion distance
		}
		if(min_dis>dis){
			a.px=distgas(engine)*1e-5;
			a.py=distgas(engine)*1e-5;
			a.pz=distgas(engine)*1e-5;
			a.mass=pp->Mgas;

			vars->effectiveOut[1][i]=a;
			i++;
		}
		collisionFlagGas.push_back(0);
	} while(i<Nof_around_gas);
}


/////////////////////////////////////////////////////////////////////
/*
	- Randomly arraying gas molecule around an ion with avoiding the
	overlapping. The velocity is picked from  the Maxwell-Boltzumann
	distribution
	- Set ion's center of mass (maybe -> 0), make pair list for initial
	step of simulation, reset the margine size.
*/
/////////////////////////////////////////////////////////////////////

void
MD::initialization_vapor(void) {
	vars->effectiveOut[2].resize(Nof_around_vapor);
	double nion=vars->effectiveIn[0][0].inAtoms.size();
	double dx,dy,dz, d,min_gv, min_iv, min_vv, gv, iv, vv;
	gv=10, iv=10, vv=10;	/*	minimum vapor-vapor, vapor-ion, vapor-gas distance */

    // Set random number generator
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-d_size*0.5,d_size*0.5);

    // Main part, generate random x, y, z positions and calculate minimum gas-gas distance.
	int i=0;
	do {
 		Molecule_out a;
		min_gv=min_iv=min_vv=10000.0;
		a.qx=r(mt), a.qy=r(mt), a.qz=r(mt);
		if(i>0){
			for(auto &b : vars->effectiveOut[2]) {
				dx = a.qx - b.qx;
				dy = a.qy - b.qy;
				dz = a.qz - b.qz;
				adjust_periodic(dx, dy, dz, d_size);
				d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<min_vv) min_vv=d; // minimum vapor-vapor distance
			}
		}
		for(auto &b : vars->effectiveIn[0][0].inAtoms) {
			dx=a.qx-b.qx;
			dy=a.qy-b.qy;
			dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz, d_size);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_iv) min_iv=d; // minimum vapor-ion distance
		}
    for(auto &b : vars->effectiveOut[1]) {
			dx=a.qx-b.qx;
			dy=a.qy-b.qy;
			dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz, d_size);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_gv) min_gv=d; // minimum gas-vapor distance
		}
		if(min_gv>gv && min_iv>iv && min_vv>vv){
			a.px=distvapor(engine)*1e-5;
			a.py=distvapor(engine)*1e-5;
			a.pz=distvapor(engine)*1e-5;
			a.mass=pp->Mvapor;

			vars->effectiveOut[2][i]=a;
			i++;
		}
		collisionFlagVapor.push_back(0);
	} while(i<Nof_around_vapor);
}
