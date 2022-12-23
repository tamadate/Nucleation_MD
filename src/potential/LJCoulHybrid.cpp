#include "potential.hpp"

/////////////////////////////////////////////////////////////////////
/* Calculate force working on ion-gas (LJ) */
/////////////////////////////////////////////////////////////////////
void
PotentialLJCoulHybrid::compute(Variables *vars, FLAG *flags) {
	Molecule *mols = vars->Molecules.data();
	vars->times.tvv-=omp_get_wtime();
	for(auto &p : vars->pairsLJCoulHybrid){
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
			double force_pair = force_lj*r2inv*(1-w);
			mols[i].fx += force_pair * dx;
			mols[i].fy += force_pair * dy;
			mols[i].fz += force_pair * dz;
			mols[j].fx -= force_pair * dx;
			mols[j].fy -= force_pair * dy;
			mols[j].fz -= force_pair * dz;
			if(flags->eflag) {
				vars->Utotal.Uvv+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0)*(1-w);
			}
		}
		if((UNION&AA)==AA){
			for (auto &av1 : mols[i].inAtoms){
				for (auto &av2 : mols[j].inAtoms){
					double dx = av2.qx - av1.qx;
					double dy = av2.qy - av1.qy;
					double dz = av2.qz - av1.qz;
					double rsq = (dx * dx + dy * dy + dz * dz);
					double r2inv = 1/rsq;
					double r6inv = r2inv * r2inv * r2inv;
					int type1=av1.type;
					int type2=av2.type;
					double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
					double force_coul = qqrd2e * av1.charge * av2.charge * sqrt(r2inv);
					double force_pair = (force_lj + force_coul)*r2inv*w;
					av2.fx += force_pair * dx;
					av2.fy += force_pair * dy;
					av2.fz += force_pair * dz;
					av1.fx -= force_pair * dx;
					av1.fy -= force_pair * dy;
					av1.fz -= force_pair * dz;
					if(flags->eflag) {
						vars->Utotal.Uvv+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0)*w;
						vars->Utotal.Uvv+=force_coul*w;
					}
				}
			}
		}
	}
	vars->times.tvv+=omp_get_wtime();
}
