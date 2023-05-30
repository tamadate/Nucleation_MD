#include "../potential.hpp"

void
PotentialGasIon::compute(Variables *vars) {
	Molecule *gases = vars->gases.data();
	Atom *ions = vars->ions.data();
	vars->times.tgi-=omp_get_wtime();
	for(auto &p : pairs){
		int i=p.i;
		int j=p.j;
		for (auto &ag : gases[i].inAtoms){
			double dx = ag.qx - ions[j].qx;
			double dy = ag.qy - ions[j].qy;
			double dz = ag.qz - ions[j].qz;
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			int type1=ions[j].type;
			int type2=ag.type;
			double r6inv = r2inv * r2inv * r2inv;
			double force_pair = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1])*r2inv;
			ag.fx += force_pair * dx;
			ag.fy += force_pair * dy;
			ag.fz += force_pair * dz;
			ions[j].fx -= force_pair * dx;
			ions[j].fy -= force_pair * dy;
			ions[j].fz -= force_pair * dz;
			if(vars->eflag) {
				vars->U.Ugi+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			}
				//	if(vars->eflag) vars->totalPotential+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
				//	vars->totalVirial+=force_lj;
		}
	}
	vars->times.tgi+=omp_get_wtime();
}

void
PotentialGasIon::makePair(Variables *vars){
	pairs.clear();
	Molecule *gases = vars->gases.data();
	Atom *ions = vars->ions.data();
	for (auto i : vars->gas_in){
		for(int j=0;j<vars->ions.size();j++){
			double dx=gases[i].qx-ions[j].qx;
			double dy=gases[i].qy-ions[j].qy;
			double dz=gases[i].qz-ions[j].qz;
			double r2 = (dx * dx + dy * dy + dz * dz);
			if(r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				pairs.push_back(p);
			}
		}
	}
}