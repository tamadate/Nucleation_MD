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
	- Randomly aranging gas molecule around an ion with avoiding the
	overlapping.
	- The velocity is picked from  the Maxwell-Boltzumann
	distribution
*/
/////////////////////////////////////////////////////////////////////

void
MD::initialization_gas(void) {
	vars->CG[1].resize(Nof_around_gas); // reserve memorry
	double dx,dy,dz, d;
	double dis=5;	/*	minimum gas-gas, gas-ion distance */

  // Maxwell-Boltzumann distribution generator
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-d_size*0.5,d_size*0.5);

  // Sampling gas molecules
	int i=0;
	do {
		Molecule a;
		double min_dis=10000.0;
		// sample random position (x, y, z)
		a.qx=r(mt), a.qy=r(mt), a.qz=r(mt);

		// calculate minimum distance from existing atoms (min_dis)
		// min distance from existing gas molecules
		if(i>0){
			for(auto &b : vars->CG[1]){
				dx=a.qx-b.qx;
				dy=a.qy-b.qx;
				dz=a.qx-b.qx;
				adjust_periodic(dx, dy, dz, d_size);
				d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<min_dis) min_dis=d; // minimum gas-gas distance
			}
		}
		// min distance from existing gas molecules
		for(auto &b : vars->CG[0][0].inAtoms) {
			dx=a.qx-b.qx;
			dy=a.qy-b.qy;
			dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz, d_size);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_dis) min_dis=d; // minimum gas-ion distance
		}

		// if min_dis is less than criteria, add the sampled molecule (a)
		if(min_dis>dis){
			// sample velocities from Maxwell-Boltzumann distribution
			// set mass and id
			a.px=distgas(engine)*1e-5;
			a.py=distgas(engine)*1e-5;
			a.pz=distgas(engine)*1e-5;
			a.mass=pp->Mgas;
			a.id=i;

			vars->CG[1][i]=a;
			i++;
		}
		//collisionFlagGas.push_back(0);
	} while(i<Nof_around_gas);
}
