#include "trajectory.hpp"

// This function records the minimum distance of binding vapor molecule with ion

void 
Trajectory::postPosition(void){
    int Nion=vars->ions.size(); // number of atoms in the cneter ion
    // loop for vapor molecule in all-atom region
	// where i is the vapor molecule id and j is the atom id in the ion 
    // this loop calculate atom pair minimum distance from atoms in the center ion for each vapor molecule
    for (auto &i : vars->vapor_in) {
        double dmin=10000000;   // minimum distance
        int jmin=-1;            // atom index showing minimum distance 

        // this loop calcualte minimum vapor-ion distance (inter-atomic)
        for (auto &b : vars->vapors[i].inAtoms){   
            for (int j=0;j<Nion;j++) {
                double dx = b.qx - vars->ions[j].qx;
                double dy = b.qy - vars->ions[j].qy;
                double dz = b.qz - vars->ions[j].qz;
                double r2 = (dx * dx + dy * dy + dz * dz);
                double r = sqrt(r2);
                // collision distance defined from L-J potential
                double sig=sigma216[b.type][vars->ions[j].type]; 
                double d=r-sig;
                // if d is smaller than minimum distance so far, 
			    // the minimum distance and the closest atom ID are overwritten 
                if(dmin>d) {dmin=d; jmin=j;}
            }
        }

        // if vapor status is "binding", the results are exported
        if (collisionFlagVapor[i]!=0) {
            FILE*f=fopen(filename, "a");
            fprintf(f, "%d\t%e\t%e\t%d\n", i,vars->time,dmin,jmin);
            fclose(f);
            // and if distance is larger than the leaving distance,
            // the vapor status is switched back to "no-binding"
            if (dmin > dl) collisionFlagVapor[i]=0;
        }
        // if vapor status is "no-binding"
        else{
            // and if distance is smaller than the arrival (collision) distance,
            // the vapor status is switched to "binding"
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