#include "../md.hpp"

void
Observer::outputIonCenter(void){
	FILE*f=fopen(fileIonCenter, "a");
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, vars->IonX[0], vars->IonX[1], vars->IonX[2], vars->IonV[0], vars->IonV[1], vars->IonV[2]);
	fclose(f);
}

void
Observer::outputGasCenter(void){
	FILE*f=fopen(fileGasCenter, "a");
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, vars->gasX[0], vars->gasX[1], vars->gasX[2], vars->gasV[0], vars->gasV[1], vars->gasV[2]);
	fclose(f);
}

