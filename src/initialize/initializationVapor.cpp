#include "../md.hpp"

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
MD::initializatIonVapor(void) {
	double nion=vars->ions.size();
	double dx,dy,dz, d,min_gv, min_iv, min_vv, gv, iv, vv;
	gv=10, iv=10, vv=10;	/*	minimum vapor-vapor, vapor-ion, vapor-gas distance */

    // Set random number generator
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distvapor(0.0, sqrt(kb*pp->T/pp->mvapor));
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-con->HL,con->HL);

    // Main part, generate random x, y, z positions and calculate minimum gas-gas distance.
	int i=0;
	do {
 		Molecule a;
		min_gv=min_iv=min_vv=10000.0;
		a.qx=r(mt), a.qy=r(mt), a.qz=r(mt);
		if(i>0){
			for(auto &b : vars->vapors) {
				dx = a.qx - b.qx;
				dy = a.qy - b.qy;
				dz = a.qz - b.qz;
				adjust_periodic(dx, dy, dz, con->L);
				d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<min_vv) min_vv=d; // minimum vapor-vapor distance
			}
		}
		for(auto &b : vars->ions) {
			dx=a.qx-b.qx;
			dy=a.qy-b.qy;
			dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz, con->L);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_iv) min_iv=d; // minimum vapor-ion distance
		}
        for(auto &b : vars->gases) {
			dx=a.qx-b.qx;
			dy=a.qy-b.qy;
			dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz, con->L);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_gv) min_gv=d; // minimum gas-vapor distance
		}
		if(min_gv>gv && min_iv>iv && min_vv>vv){
			a.px=distvapor(engine)*1e-5;
			a.py=distvapor(engine)*1e-5;
			a.pz=distvapor(engine)*1e-5;
			a.mass=pp->Mvapor;
			a.inAtoms=vars->atomVapor;
			a.bonds=vars->bonds_v;
			a.angles=vars->angles_v;
			a.dihedrals=vars->dihedrals_v;
			a.inFlag=0;
			vars->vapors.push_back(a);
			i++;
		}
	} while(i<con->Nof_around_vapor);
}
