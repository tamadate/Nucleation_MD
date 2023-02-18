#include "potential.hpp"


/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on gas-gas (LJ)
	1 Td = E/N = 10e-17 Vcm2 = 10-21 Vm2
	N = 1e5 / (1.38e-23 * 300) = 2.41e25 #/m3
	E = N * 10-21 = 2.41e4 V/m
	F = qE = 1.6e-19 * 2.41e4 * z = 3.856e-15 * z J/m
		= 3.856e-15 * 1e-3 / 4.18 * 6.02e23 * 1e-10 * z kcal/mol/A
		= 5.55e-5 * z kcal/mol/A
*/
/////////////////////////////////////////////////////////////////////
void
PotentialEfield::compute(Variables *vars, FLAG *flags) {
    for (auto &a : vars->ions) {
		a.fx+=5.55e-5*a.charge*Ecoeff[0];
		a.fy+=5.55e-5*a.charge*Ecoeff[1];
		a.fz+=5.55e-5*a.charge*Ecoeff[2];
	}

}


/////////////////////////////////////////////////////////////////////
/*
	- Calculate the ion induced dipole force working between ion and gas
	molecules. Because ion have a charge, ion induce the gas dipole.
*/
/////////////////////////////////////////////////////////////////////

// charge devided by number of atoms in ions
void
PotentialIonDipole::compute(Variables *vars, FLAG *flags) {
	Molecule *gases = vars->gases.data();
	const int gs = vars->gas_in.size();
	const int is = vars->ions.size();

    /*for (int k = 0; k < gs; k++) {
        for (auto &a : vars->ions) {
            int i = vars->gas_in[k];
            double dx = gases[i].qx - a.qx;
            double dy = gases[i].qy - a.qy;
            double dz = gases[i].qz - a.qz;
			adjust_periodic(dx, dy, dz);
			double rsq = (dx * dx + dy * dy + dz * dz);

            double r2inv = 1/rsq;
            double ion_dipole = -2*alphagas*r2inv*r2inv*qqrd2e/is/is*zion*zion;
            double force_pair=ion_dipole*r2inv;
            gases[i].fx += force_pair * dx;
            gases[i].fy += force_pair * dy;
            gases[i].fz += force_pair * dz;
            a.fx -= force_pair * dx;
            a.fy -= force_pair * dy;
            a.fz -= force_pair * dz;
			if(flags->eflag) vars->totalPotential += -0.5*alphagas*r2inv*r2inv*qqrd2e/is/is;
		}
	}*/
}

