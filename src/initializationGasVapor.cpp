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
	double gas[Nof_around_gas][3], dx,dy,dz, d,min_gg,min_gi, gasgas, gasion;
	gasgas=5, gasion=5;	/*	minimum gas-gas, gas-ion distance */
    
    // Set random number generator
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-d_size*0.5,d_size*0.5);
    
    // Main part, generate random x, y, z positions and calculate minimum gas-gas distance.
	for(int i=0;i<Nof_around_gas;i++)	gas[i][0]=gas[i][1]=gas[i][2]=0;
	int i=0;
	do {
		min_gg=min_gi=10000.0;
		gas[i][0]=r(mt), gas[i][1]=r(mt), gas[i][2]=r(mt);
		for(int j=0;j<i;j++) {
			dx=gas[i][0]-gas[j][0];
			dy=gas[i][1]-gas[j][1];
			dz=gas[i][2]-gas[j][2];
			adjust_periodic(dx, dy, dz);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_gg) min_gg=d; // minimum gas-gas distance
		}
		for(auto &a : vars->ions) {
			dx=gas[i][0]-a.qx;
			dy=gas[i][1]-a.qy;
			dz=gas[i][2]-a.qz;
			adjust_periodic(dx, dy, dz);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_gi) min_gi=d; // minimum gas-ion distance
		}
		if(min_gg>gasgas && min_gi>gasion){
			double vx=distgas(engine)*1e-5;
			double vy=distgas(engine)*1e-5;
			double vz=distgas(engine)*1e-5;
			if (gastype==1)	vars->add_gases(nion+i+1, 1, gas[i][0], gas[i][1], gas[i][2], vx, vy, vz, 0.0, 0.0, 0.0, 0.0, 4.0027);
            if (gastype==2)  vars->add_gases(nion+i+1, 2, gas[i][0], gas[i][1], gas[i][2], vx, vy, vz, 0.0, 0.0, 0.0, 0.0, 28.02);
			if (gastype==3)	vars->add_gases(nion+i+1, 3, gas[i][0], gas[i][1], gas[i][2], vx, vy, vz, 0.0, 0.0, 0.0, 0.0, 28.02);
			if (gastype==4)	vars->add_gases(nion+i+1, 4, gas[i][0], gas[i][1], gas[i][2], vx, vy, vz, 0.0, 0.0, 0.0, 0.0, 39.948);
			i++;
		}
		collisionFlagGas.push_back(0);
	} while(i<Nof_around_gas);
    if (gastype==2){
        for(i=0; i<Nof_around_gas; i++){
            vars->add_gases(nion+i+1+Nof_around_gas, 2, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0);
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
	double vapor[Nof_around_vapor][3], dx,dy,dz, d,min_gv, min_iv, min_vv, gv, iv, vv;
	gv=10, iv=10, vv=10;	/*	minimum vapor-vapor, vapor-ion, vapor-gas distance */
    
    // Set random number generator
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-d_size*0.5,d_size*0.5);
    
    // Main part, generate random x, y, z positions and calculate minimum gas-gas distance.
	for(int i=0;i<Nof_around_vapor;i++)	vapor[i][0] = vapor[i][1] = vapor[i][2] = 0;
	int i=0;
	do {
		min_gv=min_iv=min_vv=10000.0;
		vapor[i][0]=r(mt), vapor[i][1]=r(mt), vapor[i][2]=r(mt);
		for(int j=0;j<i;j++) {
			dx = vapor[i][0] - vapor[j][0];
			dy = vapor[i][1] - vapor[j][1];
			dz = vapor[i][2] - vapor[j][2];
			adjust_periodic(dx, dy, dz);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_vv) min_vv=d; // minimum gas-gas distance
		}
		for(auto &a : vars->ions) {
			dx=vapor[i][0]-a.qx;
			dy=vapor[i][1]-a.qy;
			dz=vapor[i][2]-a.qz;
			adjust_periodic(dx, dy, dz);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_iv) min_iv=d; // minimum gas-ion distance
		}
        for(auto &a : vars->gases) {
			dx=vapor[i][0]-a.qx;
			dy=vapor[i][1]-a.qy;
			dz=vapor[i][2]-a.qz;
			adjust_periodic(dx, dy, dz);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_gv) min_gv=d; // minimum gas-ion distance
		}
		if(min_gv>gv && min_iv>iv && min_vv>vv){
			double vx=distvapor(engine)*1e-5;
			double vy=distvapor(engine)*1e-5;
			double vz=distvapor(engine)*1e-5;
			vars->add_vapors(nion+Nof_around_gas+i+1, 100+vaportype, vapor[i][0], vapor[i][1], vapor[i][2], vx, vy, vz, pp->Mvapor);
			i++;
		}
		collisionFlagVapor.push_back(0);
	} while(i<Nof_around_vapor);
}


