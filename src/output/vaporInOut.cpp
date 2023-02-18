#include "../md.hpp"

void
Observer::Ovin(int i, double time){
	FILE*f=fopen(fileVaporIn, "a");
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,time,vars->vapors[i].qx-vars->ion_r[0],vars->vapors[i].qy-vars->ion_r[1],vars->vapors[i].qz-vars->ion_r[2],
    vars->vapors[i].px-vars->ion_v[0],vars->vapors[i].py-vars->ion_v[1],vars->vapors[i].pz-vars->ion_v[2]);
	fclose(f);
}

void
Observer::Ovout(int i, double time){
	FILE*f=fopen(fileVaporOut, "a");
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,time,vars->vapors[i].qx-vars->ion_r[0],vars->vapors[i].qy-vars->ion_r[1],vars->vapors[i].qz-vars->ion_r[2],
        vars->vapors[i].px-vars->ion_v[0],vars->vapors[i].py-vars->ion_v[1],vars->vapors[i].pz-vars->ion_v[2]);
	fclose(f);
}