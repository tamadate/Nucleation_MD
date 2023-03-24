#include "../md.hpp"

void
Collision::outputIn(int i){
	FILE*f=fopen(fileVaporIn, "a");
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,time,vars->vapors[i].qx-vars->IonX[0],vars->vapors[i].qy-vars->IonX[1],vars->vapors[i].qz-vars->IonX[2],
    vars->vapors[i].px-vars->IonV[0],vars->vapors[i].py-vars->IonV[1],vars->vapors[i].pz-vars->IonV[2]);
	fclose(f);
}

void
Collision::outputOut(int i){
	FILE*f=fopen(fileVaporOut, "a");
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,time,vars->vapors[i].qx-vars->IonX[0],vars->vapors[i].qy-vars->IonX[1],vars->vapors[i].qz-vars->IonX[2],
        vars->vapors[i].px-vars->IonV[0],vars->vapors[i].py-vars->IonV[1],vars->vapors[i].pz-vars->IonV[2]);
	fclose(f);
}