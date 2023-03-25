#include "../potential.hpp"


/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on ion-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialVaporIon::compute(Variables *vars) {
	Molecule *vapors = vars->vapors.data();
	Atom *ions = vars->ions.data();
	const int is = vars->ions.size();
	const int vin_size = vars->vapor_in.size();
	vars->times.tvi-=omp_get_wtime();
	for(int i=0;i<vin_size;i++){
		for (auto &av : vapors[vars->vapor_in[i]].inAtoms){
			for(auto &ai : vars->ions){
				double dx = av.qx - ai.qx;
				double dy = av.qy - ai.qy;
				double dz = av.qz - ai.qz;
				double rsq = (dx * dx + dy * dy + dz * dz);
				double r2inv = 1/rsq;
				int type1=av.type;
				int type2=ai.type;
				double r6inv = r2inv * r2inv * r2inv;
				double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
				double force_coul = qqrd2e * av.charge * ai.charge * sqrt(r2inv);
				double force_pair = (force_lj + force_coul)*r2inv;
				av.fx += force_pair * dx;
				av.fy += force_pair * dy;
				av.fz += force_pair * dz;
				ai.fx -= force_pair * dx;
				ai.fy -= force_pair * dy;
				ai.fz -= force_pair * dz;
				if(vars->eflag) {
					vars->U.Uvi+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
					vars->U.Uvi+=force_coul;
				}
				//	vars->totalVirial+=force_lj;
			}
		}
	}
	vars->times.tvi+=omp_get_wtime();
}


void
PotentialVaporIon::makePair(Variables *vars) {}