#include "stickPosition.hpp"

void
StickPosition::postLoop(void){
	for(auto &i : vars->vapor_in){
		double minDist=1e10;
		int closeAtomID=0;
		for(auto &j : stickPositionList){
			double dx=vars->vapors[i].qx-vars->IonX[0];
			double dy=vars->vapors[i].qy-vars->IonX[1];
			double dz=vars->vapors[i].qz-vars->IonX[2];
			double dr2=dx*dx+dy*dy+dz*dz;
			if(dr2<minDist){
				closeAtomID=j;
				minDist=dr2;
			}
		}
		FILE*f=fopen(filename, "a");
		fprintf(f, "%e\t%d\t%d\t%e\n", vars->time,i,closeAtomID,sqrt(minDist));
		fclose(f);
	}
}


void 
StickPosition::initial(Variables *VARS, Physical *PP, MDcondition *CON, Observer *OBS){
	vars=VARS;
	pp=PP;
	con=CON;
	obs=OBS;

	// vapor in/out files initialization
	sprintf(filename, "stick_%d.dat", int(vars->calcID));
	FILE*f=fopen(filename, "w");
	fclose(f);	
}

