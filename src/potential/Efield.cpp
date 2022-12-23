#include "potential.hpp"

/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on gas-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialEfield::compute(Variables *vars, FLAG *flags) {
    for (auto &a : vars->Molecules[0].inAtoms) {
		a.fx+=6.2665e-5*a.charge*Ecoeff[0];
		a.fy+=6.2665e-5*a.charge*Ecoeff[1];
		a.fz+=6.2665e-5*a.charge*Ecoeff[2];
	}

}
