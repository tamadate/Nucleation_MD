//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------

void
MD:: updateGasinCenters(void){
	for (auto i : vars->gas_in){
		double X=0;
		double Y=0;
		double Z=0;
		double VX=0;
		double VY=0;
		double VZ=0;
		double Mass=0;
		for (auto &a : vars->gases[i].inAtoms){
			X+=a.qx*a.mass;
			Y+=a.qy*a.mass;
			Z+=a.qz*a.mass;
			VX+=a.px*a.mass;
			VY+=a.py*a.mass;
			VZ+=a.pz*a.mass;
			Mass+=a.mass;
		}
			//  averaged position
		vars->gases[i].qx=X/Mass;
		vars->gases[i].qy=Y/Mass;
		vars->gases[i].qz=Z/Mass;
		//  averaged velocity
		vars->gases[i].px=VX/Mass;
		vars->gases[i].py=VY/Mass;
		vars->gases[i].pz=VZ/Mass;
	}
}

void
MD:: updateVaporinCenters(void){
	for (auto i : vars->vapor_in){
		double X=0;
		double Y=0;
		double Z=0;
		double VX=0;
		double VY=0;
		double VZ=0;
		double Mass=0;
		for (auto &a : vars->vapors[i].inAtoms){
			X+=a.qx*a.mass;
			Y+=a.qy*a.mass;
			Z+=a.qz*a.mass;
			VX+=a.px*a.mass;
			VY+=a.py*a.mass;
			VZ+=a.pz*a.mass;
			Mass+=a.mass;
		}
        //  averaged position
	    vars->vapors[i].qx=X/Mass;
	    vars->vapors[i].qy=Y/Mass;
	    vars->vapors[i].qz=Z/Mass;
        //  averaged velocity
        vars->vapors[i].px=VX/Mass;
        vars->vapors[i].py=VY/Mass;
        vars->vapors[i].pz=VZ/Mass;
	}
}
