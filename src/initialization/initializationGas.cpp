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
	double dx,dy,dz, d;
	double dis=5;	/*	minimum gas-gas, gas-ion distance */
	int Nsofar=vars->Molecules.size();

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
		int loop=0;
		for(auto &b : vars->Molecules){
			// min distance from ion atoms
			if(loop==0){
				for(auto &c : b.inAtoms) {
					double dx=a.qx-c.qx;
					double dy=a.qy-c.qy;
					double dz=a.qz-c.qz;
					adjust_periodic(dx, dy, dz, d_size);
					double d=sqrt(dx*dx+dy*dy+dz*dz);
					if(d<min_dis) min_dis=d; // minimum gas-ion distance
				}
			}
			// min distance from existing gas & vapor molecules
			else{
				double dx=a.qx-b.qx;
				double dy=a.qy-b.qx;
				double dz=a.qx-b.qx;
				adjust_periodic(dx, dy, dz, d_size);
				double d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<min_dis) min_dis=d; // minimum gas-gas distance
			}
			loop++;
		}

		// if min_dis is less than criteria, add the sampled molecule (a)
		if(min_dis>dis){
			// sample velocities from Maxwell-Boltzumann distribution
			// set mass and id
			a.px=distgas(engine)*1e-5;
			a.py=distgas(engine)*1e-5;
			a.pz=distgas(engine)*1e-5;
			a.mass=pp->Mgas;
			a.type=1;
			a.id=i;
			vars->Molecules.push_back(a);
			vars->MolID[1].push_back(Nsofar+i);
			i++;
		}
		//collisionFlagGas.push_back(0);
	} while(i<Nof_around_gas);
}
