#include "potential.hpp"
/*########################################################################################

Compute ion intramolecular interaction _ Born-Mayer-Huggins-Fumi-Tosi
- Currently, NaCl is only avairable particle

#######################################################################################*/

/**********************************Force calculation******************************************/
void
PotentialBorn::compute(Variables *vars, FLAG *flags) {
	Atom *ions = vars->Molecules[0].inAtoms.data();
	const int is = vars->Molecules[0].inAtoms.size();
	vars->times.tion-=omp_get_wtime();
	for(int i=0; i<is-1; i++){
		for(int j=i+1; j<is; j++){
			double dx = ions[i].qx - ions[j].qx;
			double dy = ions[i].qy - ions[j].qy;
			double dz = ions[i].qz - ions[j].qz;
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			double r6inv = r2inv * r2inv * r2inv;
			double r=sqrt(rsq);
			int type1=ions[i].type;
			int type2=ions[j].type;

//vars->bornCoeff[type1][type2][0]=A
//vars->bornCoeff[type1][type2][1]=6C
//vars->bornCoeff[type1][type2][2]=8D
//vars->bornCoeff[type1][type2][3]=sigma
//vars->bornCoeff[type1][type2][4]=1/rho

			double rexp=exp(vars->bornCoeff[type1][type2][4]*(vars->bornCoeff[type1][type2][3]-r));
			double force1 = vars->bornCoeff[type1][type2][4]*vars->bornCoeff[type1][type2][0]*r*rexp;
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
			if(flags->eflag) {
				vars->Utotal.Uion+=force_coul;
				vars->Utotal.Uion+=rexp*vars->bornCoeff[type1][type2][0];
				vars->Utotal.Uion-=vars->bornCoeff[type1][type2][1]/6.0*r6inv;
				vars->Utotal.Uion-=vars->bornCoeff[type1][type2][2]/8.0*r6inv*r2inv;
			}
			//vars->totalVirial+=force_lj;
		}
	}
	vars->times.tion+=omp_get_wtime();
}
