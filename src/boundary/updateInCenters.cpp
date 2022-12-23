#include "../md.hpp"

void
MD:: updateInCenters(void){
	int i=0;
	for (auto &a : vars->Molecules){
		if(vars->Region[i]==00000010) continue;
		double X=0;
		double Y=0;
		double Z=0;
		double VX=0;
		double VY=0;
		double VZ=0;
		double Mass=0;
		for (auto &c : a.inAtoms){
			X+=c.qx*c.mass;
			Y+=c.qy*c.mass;
			Z+=c.qz*c.mass;
			VX+=c.px*c.mass;
			VY+=c.py*c.mass;
			VZ+=c.pz*c.mass;
			Mass+=c.mass;
		}
			//  averaged position
		a.qx=X/Mass;
		a.qy=Y/Mass;
		a.qz=Z/Mass;
		//  averaged velocity
		a.px=VX/Mass;
		a.py=VY/Mass;
		a.pz=VZ/Mass;
		i++;
	}
}
