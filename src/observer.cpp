//------------------------------------------------------------------------
#include "observer.hpp"
//------------------------------------------------------------------------
void
Observer::computeProps(Variables *vars,int molID){
	Kin[molID]=0;
	double Nin=0;
	Kout[molID]=0;
	double Nout=0;

	for(auto &mg : vars->effectiveIn[molID]){
		for (auto &ag : mg.inAtoms){
			Kin[molID] += ag.px * ag.px * ag.mass;
			Kin[molID] += ag.py * ag.py * ag.mass;
			Kin[molID] += ag.pz * ag.pz * ag.mass;
			Nin++;
		}
	}
	Kin[molID]*= (0.5 * real_to_kcalmol);
	Tin[molID]=Kin[1]/double(Nin)*coeff;

	for(auto &mg : vars->effectiveOut[1]){
		Kout[molID] += mg.px * mg.px * mg.mass;
		Kout[molID] += mg.py * mg.py * mg.mass;
		Kout[molID] += mg.pz * mg.pz * mg.mass;
		Nout++;
	}
	Kout[molID]*= (0.5 * real_to_kcalmol);
	Tout[molID]=Kout[molID]/double(Nout)*coeff;

};



double
Observer::pressure(Variables *vars, std::vector<Pair> &pairs, double Treal, double virial,double p,double T) {
	double phi = 0.0;
	/*const int ps = pairs.size();
	Gas *gases = vars->gases.data();
	for (int k = 0; k < ps; k++) {
		const int i = pairs[k].i;
		const int j = pairs[k].j;
		double dx = gases[j].qx - gases[i].qx;
		double dy = gases[j].qy - gases[i].qy;
		double dz = gases[j].qz - gases[i].qz;
		adjust_periodic(dx, dy, dz);
		double r2 = (dx * dx + dy * dy + dz * dz);
		double r2inv= 1/r2;
		int type1=gases[i].type;
		int type2=gases[j].type;
		if (r2 < CL2){
			double r6inv = r2inv * r2inv * r2inv;
			phi += r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
		}
	}*/
	phi = phi * Cpress;
	return  p/T*Treal + phi;
}
