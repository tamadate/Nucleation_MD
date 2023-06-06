#include "stickPosition.hpp"

// This function records sticking position of vapors on the ion.

void
StickPosition::postLoop(void){
	// loop for vapor molecule in all-atom region
	// this calculates the distance of ion and vapor
	// where i is the vapor molecule id and j is the atom id in the ion
	for(auto &i : vars->vapor_in){
		double minDist=1e10;	// define minimum distance 
		int closeAtomID=0;		// closest atom index
		for(auto &j : stickPositionList){
			// calculate distance of i and j, dr2
			double dx=vars->vapors[i].qx-vars->ions[j].qx;
			double dy=vars->vapors[i].qy-vars->ions[j].qy;
			double dz=vars->vapors[i].qz-vars->ions[j].qz;
			double dr2=dx*dx+dy*dy+dz*dz;
			// if dr2 is smaller than minimum distance so far, 
			// the minimum distance and the closest atom ID are ovewritten 
			if(dr2<minDist){
				closeAtomID=j;
				minDist=dr2;
			}
		}
		// export values
		FILE*f=fopen(filename, "a");
		fprintf(f, "%e\t%d\t%d\t%e\n", vars->time,i,closeAtomID,sqrt(minDist));
		fclose(f);
	}
}


void 
StickPosition::initial(Variables *VARS, Physical *PP, MDcondition *CON, Observer *OBS){
	// copy pointers
	vars=VARS;
	pp=PP;
	con=CON;
	obs=OBS;

	// vapor in/out files initialization
	sprintf(filename, "stick_%d.dat", int(vars->calcID));
	FILE*f=fopen(filename, "w");
	fclose(f);	
}

