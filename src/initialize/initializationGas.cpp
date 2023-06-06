#include "../md.hpp"

/////////////////////////////////////////////////////////////////////
/*
	- Randomly aranging gas molecule around an ion with avoiding the
	overlapping. The velocity is picked from  the Maxwell-Boltzumann
	distribution
	- Set ion's center of mass (maybe -> 0), make pair list for initial
	step of simulation, reset the margine size.
*/
/////////////////////////////////////////////////////////////////////

void
MD::initialization_gas(void) {
	double nion=vars->ions.size();
	double dx,dy,dz, d;
	double dis=5;	/*	minimum gas-gas, gas-ion distance */

    // Set random number generator
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distgas(0.0, sqrt(kb*pp->T/pp->mgas));
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-con->HL,con->HL);

    // Main part, generate random x, y, z positions and calculate minimum gas-gas distance.
	int i=0;
	do {
		Molecule a;
		double min_dis=10000.0;
		a.qx=r(mt), a.qy=r(mt), a.qz=r(mt);
		if(i>0){
			for(auto &b : vars->gases){
				dx=a.qx-b.qx;
				dy=a.qy-b.qx;
				dz=a.qx-b.qx;
				adjust_periodic(dx, dy, dz, con->L);
				d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<min_dis) min_dis=d; // minimum gas-gas distance
			}
		}
		for(auto &b : vars->ions) {
			dx=a.qx-b.qx;
			dy=a.qy-b.qy;
			dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz, con->L);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_dis) min_dis=d; // minimum gas-ion distance
		}
		if(min_dis>dis){
			a.px=distgas(engine)*1e-5;
			a.py=distgas(engine)*1e-5;
			a.pz=distgas(engine)*1e-5;
			a.mass=pp->Mgas;
			a.inAtoms=vars->atomGas(gastype);

			a.inFlag=0;
			vars->gases.push_back(a);
			i++;
		}
	} while(i<con->Nof_around_gas);
}
