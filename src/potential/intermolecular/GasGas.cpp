#include "../potential.hpp"

void
PotentialGasGas::compute(Variables *vars) {
	Molecule *gases = vars->gases.data();
	for(auto &p : pairs){
		int i=p.i;
		int j=p.j;
		for (auto &ag1 : gases[i].inAtoms){
			for (auto &ag2 : gases[j].inAtoms){
				double dx = ag1.qx - ag2.qx;
				double dy = ag1.qy - ag2.qy;
				double dz = ag1.qz - ag2.qz;
				double rsq = (dx * dx + dy * dy + dz * dz);
				double r2inv = 1/rsq;
				int type1=ag1.type;
				int type2=ag2.type;
				double r6inv = r2inv * r2inv * r2inv;
				double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
				double force_coul = qqrd2e * ag1.charge * ag2.charge * sqrt(r2inv);
				double force_pair = (force_lj + force_coul)*r2inv;
				ag1.fx += force_pair * dx;
				ag1.fy += force_pair * dy;
				ag1.fz += force_pair * dz;
				ag2.fx -= force_pair * dx;
				ag2.fy -= force_pair * dy;
				ag2.fz -= force_pair * dz;
				if(vars->eflag) {
					vars->U.Ugg+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
				}
				//vars->totalVirial+=force_lj;
			}
		}
	}
}

void 
PotentialGasGas::makePair(Variables *vars){
	pairs.clear();
	Molecule *gases = vars->gases.data();
	const int gs = vars->gas_in.size();
	if(gs>1){
		Pair p;
		for(int i1=0; i1<gs-1; i1++){
			p.i=vars->gas_in[i1];
			for(int i2=i1+1; i2<gs; i2++){
				p.j=vars->gas_in[i2];
				double dx = gases[p.i].qx - gases[p.j].qx;
				double dy = gases[p.i].qy - gases[p.j].qy;
				double dz = gases[p.i].qz - gases[p.j].qz;
				double r2 = (dx * dx + dy * dy + dz * dz);
				if (r2 < ML2) pairs.push_back(p);
			}
		}
	}
}
