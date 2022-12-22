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
	Atom *ions = vars->CG[0][0].inAtoms.data();
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
		ions[i].fx += force_pair * dx;
		ions[i].fy += force_pair * dy;
		ions[i].fz += force_pair * dz;
		ions[j].fx -= force_pair * dx;
		ions[j].fy -= force_pair * dy;
		ions[j].fz -= force_pair * dz;
		if(flags->eflag) {
			vars->Utotal.Uion+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			vars->Utotal.Uion+=force_coul;
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
	for(auto &mol : vars->CG[1]){
		if(mol.inFlag==0) continue;
		Atom *g1= & mol.inAtoms[0];
		Atom *g2= & mol.inAtoms[1];
    double dx = mol.inAtoms[0].qx - mol.inAtoms[1].qx;
    double dy = mol.inAtoms[0].qy - mol.inAtoms[1].qy;
    double dz = mol.inAtoms[0].qz - mol.inAtoms[1].qz;
    double rsq = (dx * dx + dy * dy + dz * dz);
    double r = sqrt(rsq);
		double dr= r - 1.098;
    double rk = 1221.7 * dr;
    double force_bond_harmonic = -2.0*rk/r;
    mol.inAtoms[0].fx += force_bond_harmonic * dx;
    mol.inAtoms[0].fy += force_bond_harmonic * dy;
    mol.inAtoms[0].fz += force_bond_harmonic * dz;
    mol.inAtoms[1].fx -= force_bond_harmonic * dx;
    mol.inAtoms[1].fy -= force_bond_harmonic * dy;
    mol.inAtoms[1].fz -= force_bond_harmonic * dz;
		if(flags->eflag) vars->Utotal.Ugas+=rk*dr;
	}
	vars->times.tgas+=omp_get_wtime();
}
