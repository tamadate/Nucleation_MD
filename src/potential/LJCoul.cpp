//------------------------------------------------------------------------
#include "potential.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on ion-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialLJCoul::compute(Variables *vars, FLAG *flags) {
	Molecule *mols = vars->Molecules.data();
	Atom *ions = vars->Molecules[0].inAtoms.data();
	vars->times.tvi-=omp_get_wtime();
	for(auto &p : vars->pairsLJCoul){
		int i=p.i;
		int j=p.j;
		for (auto &at : mols[i].inAtoms){
			double dx = at.qx - ions[j].qx;
			double dy = at.qy - ions[j].qy;
			double dz = at.qz - ions[j].qz;
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			int type1=at.type;
			int type2=ions[j].type;
			double r6inv = r2inv * r2inv * r2inv;
			double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
			double force_coul = qqrd2e * at.charge * ions[j].charge * sqrt(r2inv);
			double force_pair = (force_lj + force_coul)*r2inv;
			at.fx += force_pair * dx;
			at.fy += force_pair * dy;
			at.fz += force_pair * dz;
			ions[j].fx -= force_pair * dx;
			ions[j].fy -= force_pair * dy;
			ions[j].fz -= force_pair * dz;
			if(flags->eflag) {
				vars->Utotal.Uvi+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
				vars->Utotal.Uvi+=force_coul;
			}
		}
	}
	vars->times.tvi+=omp_get_wtime();
}
