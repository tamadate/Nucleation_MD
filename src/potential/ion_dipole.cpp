#include "potential.hpp"

/////////////////////////////////////////////////////////////////////
/*
	- Calculate the ion induced dipole force working between ion and gas
	molecules. Because ion have a charge, ion induce the gas dipole.
*/
/////////////////////////////////////////////////////////////////////


// charge devided by number of atoms in ions
void
PotentialIonDipole::compute(Variables *vars, FLAG *flags) {
	//Molecule *gases = vars->AA[1].data();

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
