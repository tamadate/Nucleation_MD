#include "trajectory.hpp"


void 
Trajectory::postPosition(void){
    int Nion=vars->ions.size();
    for (auto &i : vars->vapor_in) {
        double dmin=10000000;
        int jmin=-1;
        for (auto &b : vars->vapors[i].inAtoms){
            for (int j=0;j<Nion;j++) {
                double dx = b.qx - vars->ions[j].qx;
                double dy = b.qy - vars->ions[j].qy;
                double dz = b.qz - vars->ions[j].qz;
                double r2 = (dx * dx + dy * dy + dz * dz);
                double r = sqrt(r2);
                double sig=sigma216[b.type][vars->ions[j].type];
                double d=r-sig;
                if(dmin>d) {dmin=d; jmin=j;}
            }
        }

        if (collisionFlagVapor[i]!=0) {
            FILE*f=fopen(filename, "a");
            fprintf(f, "%d\t%e\t%e\t%d\n", i,vars->time,dmin,jmin);
            fclose(f);
            if (dmin > dl) collisionFlagVapor[i]=0;
        }
        else{
            if (dmin < da) collisionFlagVapor[i]=1;
        }
	}
}

void 
Trajectory::initial(Variables *VARS, Physical *PP, MDcondition *CON, Observer *OBS){
	vars=VARS;
	pp=PP;
	con=CON;
	obs=OBS;
	int Nvapor=vars->vapors.size();
	collisionFlagVapor.resize(Nvapor);
	for (auto i : collisionFlagVapor) i=0;

	// vapor in/out files initialization
	sprintf(filename, "vapor_in_%d.dat", int(vars->calcID));
	FILE*f=fopen(filename, "w");
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