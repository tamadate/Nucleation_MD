//------------------------------------------------------------------------
#include "../md.hpp"
//------------------------------------------------------------------------

void
MD:: updateInCenters(void){
	for (auto &a : vars->AA){
		for (auto &b : a){
			double X=0;
			double Y=0;
			double Z=0;
			double VX=0;
			double VY=0;
			double VZ=0;
			double Mass=0;
			for (auto &c : b.inAtoms){
				X+=c.qx*c.mass;
				Y+=c.qy*c.mass;
				Z+=c.qz*c.mass;
				VX+=c.px*c.mass;
				VY+=c.py*c.mass;
				VZ+=c.pz*c.mass;
				Mass+=c.mass;
			}
				//  averaged position
			b.qx=X/Mass;
			b.qy=Y/Mass;
			b.qz=Z/Mass;
			//  averaged velocity
			b.px=VX/Mass;
			b.py=VY/Mass;
			b.pz=VZ/Mass;
		}
	}
}
