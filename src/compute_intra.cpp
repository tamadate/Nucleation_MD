//------------------------------------------------------------------------
#include "potential.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on ion-ion (Coulombic+LJ)
*/
/////////////////////////////////////////////////////////////////////
void
Potential::compute(Variables *vars, FLAG *flags){
	Atom *ions = vars->ions.data();
	vars->times.tion-=omp_get_wtime();
	int ipsize=vars->ion_pairs.size();
	#pragma omp parallel for
	for (int ip=0;ip<ipsize;ip++){
		int nth=omp_get_thread_num();
		int i=vars->ion_pairs[ip].i;
		int j=vars->ion_pairs[ip].j;
		double dx = ions[i].qx - ions[j].qx;
		double dy = ions[i].qy - ions[j].qy;
		double dz = ions[i].qz - ions[j].qz;
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r2inv = 1/rsq;
		int type1=ions[i].type;
		int type2=ions[j].type;
		double r6inv = r2inv * r2inv * r2inv;
		double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
		double force_coul = qqrd2e * ions[i].charge * ions[j].charge * sqrt(r2inv);
		double force_pair = (force_lj + force_coul)*r2inv;
		ions[i].fxMP[nth] += force_pair * dx;
		ions[i].fyMP[nth] += force_pair * dy;
		ions[i].fzMP[nth] += force_pair * dz;
		ions[j].fxMP[nth] -= force_pair * dx;
		ions[j].fyMP[nth] -= force_pair * dy;
		ions[j].fzMP[nth] -= force_pair * dz;
		if(flags->eflag) {
			vars->U_MP[nth].Uion+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			vars->U_MP[nth].Uion+=force_coul;
		}
/*		if(flags->vflag){
			vars->totalVirial+=force_lj+force_coul;
		}*/
	}
	vars->times.tion+=omp_get_wtime();
}


void
PotentialGasIntra::compute(Variables *vars, FLAG *flags){
	vars->times.tgas-=omp_get_wtime();
	for(auto &i : vars->gas_in){
		Atom *g1= & vars->gases[i].inAtoms[0];
		Atom *g2= & vars->gases[i].inAtoms[1];
    double dx = g1->qx - g2->qx;
    double dy = g1->qy - g2->qy;
    double dz = g1->qz - g2->qz;
    double rsq = (dx * dx + dy * dy + dz * dz);
    double r = sqrt(rsq);
		double dr= r - 1.098;
    double rk = 1221.7 * dr;
    double force_bond_harmonic = -2.0*rk/r;
    g1->fx += force_bond_harmonic * dx;
    g1->fy += force_bond_harmonic * dy;
    g1->fz += force_bond_harmonic * dz;
    g2->fx -= force_bond_harmonic * dx;
    g2->fy -= force_bond_harmonic * dy;
    g2->fz -= force_bond_harmonic * dz;
		if(flags->eflag) vars->Utotal.Ugas+=rk*dr;
	}
	vars->times.tgas+=omp_get_wtime();
}
