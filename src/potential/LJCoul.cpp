//------------------------------------------------------------------------
#include "potential.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on ion-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialLJ::compute(Variables *vars, FLAG *flags) {
	Molecule *mols = vars->Molecules.data();
	Atom *ions = vars->Molecules[0].inAtoms.data();
	vars->times.tgi-=omp_get_wtime();
	for(auto &p : vars->pairsLJ){
		int i=p.i;
		int j=p.j;
		for (auto &ag : mols[i].inAtoms){
			double dx = ag.qx - ions[j].qx;
			double dy = ag.qy - ions[j].qy;
			double dz = ag.qz - ions[j].qz;
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			int type1=ag.type;
			int type2=ions[j].type;
			double r6inv = r2inv * r2inv * r2inv;
			double force_pair = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1])*r2inv;
			ag.fx += force_pair * dx;
			ag.fy += force_pair * dy;
			ag.fz += force_pair * dz;
			ions[j].fx -= force_pair * dx;
			ions[j].fy -= force_pair * dy;
			ions[j].fz -= force_pair * dz;
			if(flags->eflag) {
				vars->Utotal.Ugi+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			}
				//	if(flags->eflag) vars->totalPotential+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
				//	vars->totalVirial+=force_lj;
		}
	}
	vars->times.tgi+=omp_get_wtime();
}
