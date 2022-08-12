//------------------------------------------------------------------------
#include "observer.hpp"
//------------------------------------------------------------------------
double
Observer::gas_kinetic_energy(Variables *vars) {
	double k = 0;
	for (auto &a : vars->gases) {
		k += a.px * a.px * a.mass;
		k += a.py * a.py * a.mass;
		k += a.pz * a.pz * a.mass;
	}
	return k *0.5 * real_to_kcalmol;
};

double 
Observer::gas_temperature(Variables *vars) {
	return gas_kinetic_energy(vars) / 1.5 * kb_real_inv/ static_cast<double>(Nof_around_gas);
}

double
Observer::gas_total_kinetic_energy(Variables *vars) {
	double k = 0;
	for (auto &a : vars->gases) {
		k += a.px * a.px * a.mass;
		k += a.py * a.py * a.mass;
		k += a.pz * a.pz * a.mass;
	}
	return k *0.5 * real_to_kcalmol;
};
double
Observer::ion_kinetic_energy(Variables *vars) {
	double k = 0;
	
   for (auto &a : vars->ions) {
		k += a.px * a.px * a.mass;
		k += a.py * a.py * a.mass;
		k += a.pz * a.pz * a.mass;
	}
	return k *0.5 * real_to_kcalmol;
};

double
Observer::ion_temperature(Variables *vars) {
	const int pn=vars->ions.size();
	return ion_kinetic_energy(vars) / 1.5 * kb_real_inv/ static_cast<double>(pn);
};

double
Observer::pressure(Variables *vars, std::vector<Pair> &pairs, double Treal, double virial) {
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
