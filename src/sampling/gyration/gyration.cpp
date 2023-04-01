#include "gyration.hpp"

void
Gyration::postLoop(void){
	double numerator=0;
	for(auto &a : vars->ions){
		double dx=a.qx-vars->IonX[0];
		double dy=a.qy-vars->IonX[1];
		double dz=a.qz-vars->IonX[2];
		double rsq=dx*dx+dy*dy+dz*dz;
		numerator+=a.mass*rsq;
	}
	FILE*f=fopen(filename, "a");
	fprintf(f, "%f\t%f\t",vars->time,sqrt(numerator/pp->Mion));
	fclose(f);
}


void 
Gyration::initial(Variables *VARS, Physical *PP, MDcondition *CON, Observer *OBS){
	vars=VARS;
	pp=PP;
	con=CON;
	obs=OBS;

	// vapor in/out files initialization
	sprintf(filename, "gyration_%d.dat", int(vars->calcID));
	FILE*f=fopen(filename, "w");
	fclose(f);	
}

