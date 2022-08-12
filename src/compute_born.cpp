#include "potential.hpp"
/*########################################################################################

-----compute intramolecular interaction-----

#######################################################################################*/

/**********************************Force calculation******************************************/
void
PotentialBorn::compute(Variables *vars, FLAG *flags) {
	Ion *ions = vars->ions.data();
	const int is = vars->ions.size();
	for(int i=0; i<is-1; i++){
		for(int j=i+1; j<is; j++){
			double dx = ions[i].qx - ions[j].qx;
			double dy = ions[i].qy - ions[j].qy;
			double dz = ions[i].qz - ions[j].qz;
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			double r6inv = r2inv * r2inv * r2inv;
			double r=sqrt(rsq);
			int type1=ions[i].type-11;
			int type2=ions[j].type-11;
			double force1 = vars->bornCoeff[type1][type2][0]*vars->bornCoeff[type1][type2][4]*exp(vars->bornCoeff[type1][type2][4]*(vars->bornCoeff[type1][type2][3]-r))*r;
			double force2 = -vars->bornCoeff[type1][type2][1]*r6inv;
			double force3 = -vars->bornCoeff[type1][type2][2]*r6inv*r2inv;
			double force_coul = qqrd2e*ions[i].charge*ions[j].charge*sqrt(r2inv);
			double force_pair = (force1+force2+force3+force_coul)*r2inv;

			ions[i].fx += force_pair * dx;
			ions[i].fy += force_pair * dy;
			ions[i].fz += force_pair * dz;
			ions[j].fx -= force_pair * dx;
			ions[j].fy -= force_pair * dy;
			ions[j].fz -= force_pair * dz;
			//if(flags->eflag) vars->totalPotential+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			//vars->totalVirial+=force_lj;
		}
	}
}

