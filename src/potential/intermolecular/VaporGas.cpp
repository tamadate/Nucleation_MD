#include "../potential.hpp"


/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on vapor-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialVaporGas::compute(Variables *vars, FLAG *flags) {
	Molecule *gases = vars->gases.data();
	Molecule *vapors = vars->vapors.data();
	vars->times.tvg-=omp_get_wtime();
	for(auto &p : pairs){
		int i=p.i;
		int j=p.j;
		for (auto &ag : gases[i].inAtoms){
			for (auto &av : vapors[j].inAtoms){
				double dx = av.qx - ag.qx;
				double dy = av.qy - ag.qy;
				double dz = av.qz - ag.qz;
				double rsq = (dx * dx + dy * dy + dz * dz);
				double r2inv = 1/rsq;
				int type1=ag.type;
				int type2=av.type;
				double r6inv = r2inv * r2inv * r2inv;
				double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
				double force_pair = (force_lj)*r2inv;
				av.fx += force_pair * dx;
				av.fy += force_pair * dy;
				av.fz += force_pair * dz;
				ag.fx -= force_pair * dx;
				ag.fy -= force_pair * dy;
				ag.fz -= force_pair * dz;
				if(flags->eflag) {
					vars->Utotal.Uvg+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
				}
					//	if(flags->eflag) vars->totalPotential+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
					//	vars->totalVirial+=force_lj;
			}
		}
	}
	vars->times.tvg+=omp_get_wtime();

}


void
PotentialVaporGas::makePair(Variables *vars) {
	pairs.clear();
	Molecule *vapors = vars->vapors.data();
	Molecule *gases = vars->gases.data();
	for (auto i : vars->gas_in){
		for (auto j : vars->vapor_in){
			double dx = gases[i].qx - vapors[j].qx;
			double dy = gases[i].qy - vapors[j].qy;
			double dz = gases[i].qz - vapors[j].qz;
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				pairs.push_back(p);
			}
		}
	}
}