#include "collision.hpp"


void 
Collision::postPosition(void){
	for (auto &i : vars->vapor_in) {
		double dmin=10000000;
		for (auto &b : vars->vapors[i].inAtoms){
			for (auto &a : vars->ions) {
				double dx = b.qx - a.qx;
				double dy = b.qy - a.qy;
				double dz = b.qz - a.qz;
				double r2 = (dx * dx + dy * dy + dz * dz);
				double r = sqrt(r2);
				double d=r-sigma216[b.type][a.type];
				if(dmin>d) dmin=d;
			}
		}
		if (dmin > dl && collisionFlagVapor[i]!=0) {
			FILE*f=fopen(fileVaporOut, "a");
			fprintf(f, "%d %e %e %e %e %e %e %e\n", i,vars->time,vars->vapors[i].qx-vars->IonX[0],vars->vapors[i].qy-vars->IonX[1],vars->vapors[i].qz-vars->IonX[2],
        		vars->vapors[i].px-vars->IonV[0],vars->vapors[i].py-vars->IonV[1],vars->vapors[i].pz-vars->IonV[2]);
			fclose(f);
			collisionFlagVapor[i]=0;
		}
		if (dmin < da && collisionFlagVapor[i]==0){
			FILE*f=fopen(fileVaporIn, "a");
			fprintf(f, "%d %e %e %e %e %e %e %e\n", i,vars->time,vars->vapors[i].qx-vars->IonX[0],vars->vapors[i].qy-vars->IonX[1],vars->vapors[i].qz-vars->IonX[2],
				vars->vapors[i].px-vars->IonV[0],vars->vapors[i].py-vars->IonV[1],vars->vapors[i].pz-vars->IonV[2]);
			fclose(f);
			collisionFlagVapor[i]=1;
		}
	}
}

void 
Collision::initial(Variables *VARS, Physical *PP, MDcondition *CON, Observer *OBS){
	vars=VARS;
	pp=PP;
	con=CON;
	obs=OBS;
	int Nvapor=vars->vapors.size();
	collisionFlagVapor.resize(Nvapor);
	for (auto i : collisionFlagVapor) i=0;

	// vapor in/out files initialization
	sprintf(fileVaporIn, "vapor_in_%d.dat", int(vars->calcID));
	FILE*f=fopen(fileVaporIn, "w");
	fclose(f);
	sprintf(fileVaporOut, "vapor_out_%d.dat", int(vars->calcID));
	f=fopen(fileVaporOut, "w");
	fclose(f);

	// compute collision distance 2^(1/6)*sigma
	int num_atoms=vars->atypes.size();
	sigma216.resize(num_atoms);
	for (int i=0;i<num_atoms;i++) sigma216[i].resize(num_atoms);
	for (int i=0;i<num_atoms;i++){
		for (int j=0;j<num_atoms;j++){
			double sigmai=vars->atypes[i].coeff2;
			double sigmaj=vars->atypes[j].coeff2;
			double sigmaij=(sigmai+sigmaj)*0.5;
			sigma216[i][j]=sigmaij*pow(2,1/6.0);
		}
	}
	
}
