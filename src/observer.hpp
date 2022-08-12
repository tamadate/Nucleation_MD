#pragma once
#include "variables.hpp"
//------------------------------------------------------------------------
class Observer {
public:
	double gas_kinetic_energy(Variables *vars);
	double gas_total_kinetic_energy(Variables *vars);
	double ion_kinetic_energy(Variables *vars);
	double potential_energy(Variables *vars) {return vars->totalPotential;}
	double gas_temperature(Variables *vars) ;
	double gas_total_temperature(Variables *vars) {
		return gas_total_kinetic_energy(vars) / 1.5 * kb_real_inv/double(Nof_around_gas);
	}
	double ion_temperature(Variables *vars);
	double pressure(Variables *vars, std::vector<Pair> &pairs, double Treal, double virial);
	double total_energy(Variables *vars, std::vector<Pair> &pairs) {return gas_kinetic_energy(vars) + potential_energy(vars);}
};
//------------------------------------------------------------------------
