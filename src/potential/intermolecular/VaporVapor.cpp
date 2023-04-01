#include "../potential.hpp"


/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on ion-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialVaporVapor::compute(Variables *vars) {
	Molecule *vapors = vars->vapors.data();
	const int vs = pairs.size();
	vars->times.tvv-=omp_get_wtime();
	for(int i=0;i<vs;i++){
		Pair p=pairs[i];
		for (auto &av1 : vapors[p.i].inAtoms){
			for (auto &av2 : vapors[p.j].inAtoms){
				double dx = av1.qx - av2.qx;
				double dy = av1.qy - av2.qy;
				double dz = av1.qz - av2.qz;
				double rsq = (dx * dx + dy * dy + dz * dz);
				double r2inv = 1/rsq;
				int type1=av1.type;
				int type2=av2.type;
				double r6inv = r2inv * r2inv * r2inv;
				double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
				double force_coul = qqrd2e * av1.charge * av2.charge * sqrt(r2inv);
				double force_pair = (force_lj + force_coul)*r2inv;
				av1.fx += force_pair * dx;
				av1.fy += force_pair * dy;
				av1.fz += force_pair * dz;
				av2.fx -= force_pair * dx;
				av2.fy -= force_pair * dy;
				av2.fz -= force_pair * dz;
				if(vars->eflag) {
					vars->U.Uvv+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
					vars->U.Uvv+=force_coul;
				}
				//	vars->totalVirial+=force_lj;
			}
		}
	}
	vars->times.tvv+=omp_get_wtime();
}

void
PotentialVaporVapor::makePair(Variables *vars) {
	pairs.clear();
	Molecule *vapors = vars->vapors.data();
	const int vs = vars->vapor_in.size();
	if(vs>1){
		Pair p;
		for(int i1=0; i1<vs-1; i1++){
			p.i=vars->vapor_in[i1];
			for(int i2=i1+1; i2<vs; i2++){
				p.j=vars->vapor_in[i2];
				pairs.push_back(p);
			}
		}
	}
}