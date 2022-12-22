#include "../md.hpp"


/////////////////////////////////////////////////////////////////////
/*
	- Randomly arraying gas molecule around an ion with avoiding the
	overlapping.
  - The velocity is picked from  the Maxwell-Boltzumann
	distribution
*/
/////////////////////////////////////////////////////////////////////

void
MD::initialization_vapor(void) {
	vars->CG[2].resize(Nof_around_vapor);
	double dis=10;	/*	minimum vapor-vapor, vapor-ion, vapor-gas distance */

  // Maxwell-Boltzumann distribution generator
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-d_size*0.5,d_size*0.5);

  // Sampling vapor molecules
	int i=0;
	do {
 		Molecule a;
		double minDis=10000.0;
    // sample random position (x, y, z)
		a.qx=r(mt), a.qy=r(mt), a.qz=r(mt);

    // calculate minimum distance from existing atoms (min_dis)
    // min distance from existing vapor molecules
		if(i>0){
			for(auto &b : vars->CG[2]) {
				double dx = a.qx - b.qx;
				double dy = a.qy - b.qy;
				double dz = a.qz - b.qz;
				adjust_periodic(dx, dy, dz, d_size);
				double d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<minDis) minDis=d; // minimum vapor-vapor distance
			}
		}
    // min distance from existing ion atoms
		for(auto &mol : vars->CG[0]) {
      for(auto &at : mol.inAtoms) {
  			double dx=a.qx-at.qx;
  			double dy=a.qy-at.qy;
  			double dz=a.qz-at.qz;
  			adjust_periodic(dx, dy, dz, d_size);
  			double d=sqrt(dx*dx+dy*dy+dz*dz);
  			if(d<minDis) minDis=d; // minimum vapor-ion distance
      }
		}
    // min distance from existing gas molecules
    for(auto &b : vars->CG[1]) {
			double dx=a.qx-b.qx;
			double dy=a.qy-b.qy;
			double dz=a.qz-b.qz;
			adjust_periodic(dx, dy, dz, d_size);
			double d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<minDis) minDis=d; // minimum gas-vapor distance
		}

    // if min_dis is less than criteria, add the sampled molecule (a)
		if(minDis>dis){
			a.px=distvapor(engine)*1e-5;
			a.py=distvapor(engine)*1e-5;
			a.pz=distvapor(engine)*1e-5;
			a.mass=pp->Mvapor;
      a.id=i;
      a.bonds=vars->bonds_v;
      a.angles=vars->angles_v;
      a.dihedrals=vars->dihedrals_v;
			vars->CG[2][i]=a;
			i++;
		}
		//collisionFlagVapor.push_back(0);
	} while(i<Nof_around_vapor);
}
