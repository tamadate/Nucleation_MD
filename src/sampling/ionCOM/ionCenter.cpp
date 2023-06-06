#include "ionCenter.hpp"

// Export ion center properties both of the position and velocity

void
IonCenter::postLoop(void){
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

