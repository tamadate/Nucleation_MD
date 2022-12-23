#include "potential.hpp"

/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on vapor-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialLJHybrid::compute(Variables *vars, FLAG *flags) {
	Molecule *mols = vars->Molecules.data();
	vars->times.tvg-=omp_get_wtime();
	for(auto &p : vars->pairsLJHybrid){
		int i=p.i;
		int j=p.j;
		int UNION=vars->Region[i]&vars->Region[j];	// Region[i] union Region[j]
		double w=mols[i].w*mols[j].w;
		if((UNION&CG)==CG){
			double dx = mols[i].qx - mols[j].qx;
			double dy = mols[i].qy - mols[j].qy;
			double dz = mols[i].qz - mols[j].qz;
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			double r6inv = r2inv * r2inv * r2inv;
			int type1=mols[i].type;
			int type2=mols[j].type;
			double force_lj = r6inv * (vars->pair_coeff_CG[type1][type2][0] * r6inv - vars->pair_coeff_CG[type1][type2][1]);
			double force_pair = (force_lj)*r2inv*(1-w);
			mols[i].fx += force_pair * dx;
			mols[i].fy += force_pair * dy;
			mols[i].fz += force_pair * dz;
			mols[j].fx -= force_pair * dx;
			mols[j].fy -= force_pair * dy;
			mols[j].fz -= force_pair * dz;
			if(flags->eflag) {
				vars->Utotal.Uvg+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0)*(1-w);
			}
		}
		if((UNION&AA)==AA){
			for (auto &a1 : mols[i].inAtoms){
				for (auto &a2 : mols[j].inAtoms){
					double dx = a2.qx - a1.qx;
					double dy = a2.qy - a1.qy;
					double dz = a2.qz - a1.qz;
					double rsq = (dx * dx + dy * dy + dz * dz);
					double r2inv = 1/rsq;
					double r6inv = r2inv * r2inv * r2inv;
					int type1=a1.type;
					int type2=a2.type;
					double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
					double force_pair = (force_lj)*r2inv*w;
					a2.fx += force_pair * dx;
					a2.fy += force_pair * dy;
					a2.fz += force_pair * dz;
					a1.fx -= force_pair * dx;
					a1.fy -= force_pair * dy;
					a1.fz -= force_pair * dz;
					if(flags->eflag) {
						vars->Utotal.Uvg+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0)*w;
					}
				}
			}
		}
	}
	vars->times.tvg+=omp_get_wtime();

}
