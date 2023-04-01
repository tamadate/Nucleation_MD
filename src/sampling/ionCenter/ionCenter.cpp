#include "ionCenter.hpp"

void
IonCenter::postLoop(void){
	double numerator=0;
	for(auto &a : vars->ions){
		double dx=a.qx-vars->IonX[0];
		double dy=a.qy-vars->IonX[1];
		double dz=a.qz-vars->IonX[2];
		double rsq=dx*dx+dy*dy+dz*dz;
		numerator+=a.mass*rsq;
	}
	FILE*f=fopen(filename, "a");
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, vars->IonX[0], vars->IonX[1], vars->IonX[2], vars->IonV[0], vars->IonV[1], vars->IonV[2]);
	fclose(f);
}


void 
IonCenter::initial(Variables *VARS, Physical *PP, MDcondition *CON, Observer *OBS){
	vars=VARS;
	pp=PP;
	con=CON;
	obs=OBS;

	// vapor in/out files initialization
	sprintf(filename, "ion_%d.dat", int(vars->calcID));
	FILE*f=fopen(filename, "w");
	fclose(f);	
}

