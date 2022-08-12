#include "md.hpp"
#include <random>
#include <algorithm>
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
	double nion=vars->ions.size();
	double dx,dy,dz, d,min_gg,min_gi, gasgas, gasion;
	gasgas=5, gasion=5;	/*	minimum gas-gas, gas-ion distance */
    
    // Set random number generator
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-d_size*0.5,d_size*0.5);
    
    // Main part, generate random x, y, z positions and calculate minimum gas-gas distance.
	int i=0;
	do {
		Atom a;
		a.id=nion+i+1;
		min_gg=min_gi=10000.0;
		a.qx=r(mt), a.qy=r(mt), a.qz=r(mt);
		if(i>0){
			for(auto &b : vars->gases){
				dx=a.qx-b.qx;
				dy=a.qy-b.qx;
				dz=a.qx-b.qx;
				adjust_periodic(dx, dy, dz);
				d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<min_gg) min_gg=d; // minimum gas-gas distance
			}
		}
		for(auto &b : vars->ions) {
			dx=a.qx-b.qx;
			dy=a.qy-b.qy;
			dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_gi) min_gi=d; // minimum gas-ion distance
		}
		if(min_gg>gasgas && min_gi>gasion){
			a.px=distgas(engine)*1e-5;
			a.py=distgas(engine)*1e-5;
			a.pz=distgas(engine)*1e-5;
			a.fx=a.fy=a.fz=a.charge=a.ix=a.iy=a.iz=0;
			a.type=gastype;
			a.mass=pp->Mgas;
			vars->gases.push_back(a);
			i++;
		}
		collisionFlagGas.push_back(0);
	} while(i<Nof_around_gas);
    if (gastype==2){
        for(i=0; i<Nof_around_gas; i++){
			Atom a;
			a.id=nion+i+1+Nof_around_gas;
			a.type=2;
			a.px=a.py=a.pz=a.qx=a.qy=a.qz=a.fx=a.fy=a.fz=a.charge=a.ix=a.iy=a.iz=a.mass=0;
            vars->gases.push_back(a);
        }
    }
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
	double nion=vars->ions.size();
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
 		Molecule a;
		min_gv=min_iv=min_vv=10000.0;
		a.qx=r(mt), a.qy=r(mt), a.qz=r(mt);
		if(i>0){
			for(auto &b : vars->vapors) {
				dx = a.qx - b.qx;
				dy = a.qy - b.qy;
				dz = a.qz - b.qz;
				adjust_periodic(dx, dy, dz);
				d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<min_vv) min_vv=d; // minimum gas-gas distance
			}
		}
		for(auto &b : vars->ions) {
			dx=a.qx-b.qx;
			dy=a.qy-b.qy;
			dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_iv) min_iv=d; // minimum gas-ion distance
		}
        for(auto &b : vars->gases) {
			dx=a.qx-b.qx;
			dy=a.qy-b.qy;
			dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_gv) min_gv=d; // minimum gas-ion distance
		}
		if(min_gv>gv && min_iv>iv && min_vv>vv){
			a.px=distvapor(engine)*1e-5;
			a.py=distvapor(engine)*1e-5;
			a.pz=distvapor(engine)*1e-5;
			a.mass=pp->Mvapor;
			Atom b;
			a.inAtoms=vars->atomVapor;
			a.bonds=vars->bondMeOH();
			a.angles=vars->angleMeOH();
			a.dihedrals=vars->dihedralMeOH();
			a.inFlag=0;
			vars->vapors.push_back(a);
			i++;
		}
		collisionFlagVapor.push_back(0);
	} while(i<Nof_around_vapor);
}


