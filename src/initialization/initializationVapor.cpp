#include "../md.hpp"


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
	vars->CG[2].resize(Nof_around_vapor);
	double nion=vars->AA[0][0].inAtoms.size();
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
			for(auto &b : vars->CG[2]) {
				dx = a.qx - b.qx;
				dy = a.qy - b.qy;
				dz = a.qz - b.qz;
				adjust_periodic(dx, dy, dz, d_size);
				d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<min_vv) min_vv=d; // minimum vapor-vapor distance
			}
		}
		for(auto &b : vars->AA[0][0].inAtoms) {
			dx=a.qx-b.qx;
			dy=a.qy-b.qy;
			dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz, d_size);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_iv) min_iv=d; // minimum vapor-ion distance
		}
    for(auto &b : vars->CG[1]) {
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

			vars->CG[2][i]=a;
			i++;
		}
		collisionFlagVapor.push_back(0);
	} while(i<Nof_around_vapor);
}
