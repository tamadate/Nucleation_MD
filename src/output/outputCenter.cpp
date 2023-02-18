#include "../md.hpp"

void
Observer::outputIonCenter(void){
	FILE*f=fopen(fileIonCenter, "a");
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, vars->ion_r[0], vars->ion_r[1], vars->ion_r[2], vars->ion_v[0], vars->ion_v[1], vars->ion_v[2]);
	fclose(f);
}

void
Observer::outputGasCenter(void){
	FILE*f=fopen(fileGasCenter, "a");
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, vars->gas_r[0], vars->gas_r[1], vars->gas_r[2], vars->gas_v[0], vars->gas_v[1], vars->gas_v[2]);
	fclose(f);
}

