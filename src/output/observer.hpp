#pragma once
#include "../variables.hpp"
//------------------------------------------------------------------------
class Observer {
public:
	const double coeff=kb_real_inv/1.5;
	double Kin[3];
	double Kout[3];
	double Tin[3];
	double Tout[3];


	void computeProps(Variables *vars, int molID);
	double potential_energy(Variables *vars) {return vars->totalPotential;}
	double pressure(Variables *vars, std::vector<Pair> &pairs, double Treal, double virial,double p,double T);
	//double total_energy(Variables *vars, std::vector<Pair> &pairs) {return K_g + potential_energy(vars);}
};
//------------------------------------------------------------------------
