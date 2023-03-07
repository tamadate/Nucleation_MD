//------------------------------------------------------------------------
#include "observer.hpp"
//------------------------------------------------------------------------
void
Observer::computeGasProps(Variables *vars){
	Kin_g=0;
	double Nin=0;
	Kout_g=0;
	double Nout=0;

	Molecule *gases = vars->gases.data();
	for(auto &I : vars->gas_in){
		for (auto &ag : gases[I].inAtoms){
			Kin_g += ag.px * ag.px * ag.mass;
			Kin_g += ag.py * ag.py * ag.mass;
			Kin_g += ag.pz * ag.pz * ag.mass;
			Nin++;
		}
	}
	Kin_g*= (0.5 * real_to_kcalmol);
	Tin_g=Kin_g/double(Nin)*coeff;

	for(auto i : vars->gas_out){
		Kout_g += gases[i].px * gases[i].px * gases[i].mass;
		Kout_g += gases[i].py * gases[i].py * gases[i].mass;
		Kout_g += gases[i].pz * gases[i].pz * gases[i].mass;
		Nout++;

	}
	Kout_g*= (0.5 * real_to_kcalmol);
	Tout_g=Kout_g/double(Nout)*coeff;

	K_g=Kin_g+Kout_g;
	Ngas=Nin+Nout;
	T_g=K_g/double(Ngas)*coeff;
};

void
Observer::computeIonProps(Variables *vars){
	Kion=0;
	Nion=0;
	for (auto &a : vars->ions) {
		Kion += a.px * a.px * a.mass;
		Kion += a.py * a.py * a.mass;
		Kion += a.pz * a.pz * a.mass;
		Nion++;
	}
	Kion*= (0.5 * real_to_kcalmol);
	Tion=Kion/double(Nion)*coeff;
};

void
Observer::computeVaporProps(Variables *vars){
	Kin_v=0;
	double Nin=0;
	Kout_v=0;
	double Nout=0;

	Molecule *vapors = vars->vapors.data();
	vars->vapor_in.size();
	for(auto &I : vars->vapor_in){
		for (auto &av : vapors[I].inAtoms){
			Kin_v += av.px * av.px * av.mass;
			Kin_v += av.py * av.py * av.mass;
			Kin_v += av.pz * av.pz * av.mass;
			Nin++;
		}
	}
	Kin_v*= (0.5 * real_to_kcalmol);
	if(Nin>0)	Tin_v=Kin_v/double(Nin)*coeff;
	else Tin_v=0;

	for(auto i : vars->vapor_out){
		Kout_v += vapors[i].px * vapors[i].px * vapors[i].mass;
		Kout_v += vapors[i].py * vapors[i].py * vapors[i].mass;
		Kout_v += vapors[i].pz * vapors[i].pz * vapors[i].mass;
		Nout++;
	}
	Kout_v*= (0.5 * real_to_kcalmol);
	Tout_v=Kout_v/double(Nout)*coeff;

	K_v=Kin_v+Kout_v;
	Nvap=Nin+Nout;
	T_v=K_v/double(Nvap)*coeff;
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
