#include "md.hpp"
/*########################################################################################

-----compute intramolecular interaction-----

#######################################################################################*/

/**********************************Force calculation******************************************/
void
PotentialIntraTIP3P::compute(Variables *vars, FLAG *flags) {
	Molecule *vapors = vars->vapors.data();
	Bond_type *btypes = vars->btypes.data();
	Angle_type *ctypes = vars->ctypes.data();

	double dx1, dy1, dz1, dx2, dy2, dz2, rsq1, rsq2, r1, r2, C, Cs, dtheta, tk, a, a11, a12, a22, f1[3], f3[3];

	for(auto &I : vars->vapor_in){
		dx1 = vapors[I].inAtoms[1].qx - vapors[I].inAtoms[0].qx;
		dy1 = vapors[I].inAtoms[1].qy - vapors[I].inAtoms[0].qy;
		dz1 = vapors[I].inAtoms[1].qz - vapors[I].inAtoms[0].qz;
		rsq1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
		r1 = sqrt(rsq1);
		dx2 = vapors[I].inAtoms[2].qx - vapors[I].inAtoms[0].qx;
		dy2 = vapors[I].inAtoms[2].qy - vapors[I].inAtoms[0].qy;
		dz2 = vapors[I].inAtoms[2].qz - vapors[I].inAtoms[0].qz;
		rsq2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
		r2 = sqrt(rsq2);

		double dr = (r1-0.9572);
		double rk = 450 * dr;
		double force_bond_harmonic = -rk*2/r1;
		vapors[I].inAtoms[1].fx += force_bond_harmonic * dx1;
		vapors[I].inAtoms[1].fy += force_bond_harmonic * dy1;
		vapors[I].inAtoms[1].fz += force_bond_harmonic * dz1;
		vapors[I].inAtoms[0].fx -= force_bond_harmonic * dx1;
		vapors[I].inAtoms[0].fy -= force_bond_harmonic * dy1;
		vapors[I].inAtoms[0].fz -= force_bond_harmonic * dz1;

		dr = (r2-0.9572);
		rk = 450 * dr;
		force_bond_harmonic = -rk*2/r2;
		vapors[I].inAtoms[2].fx += force_bond_harmonic * dx2;
		vapors[I].inAtoms[2].fy += force_bond_harmonic * dy2;
		vapors[I].inAtoms[2].fz += force_bond_harmonic * dz2;
		vapors[I].inAtoms[0].fx -= force_bond_harmonic * dx2;
		vapors[I].inAtoms[0].fy -= force_bond_harmonic * dy2;
		vapors[I].inAtoms[0].fz -= force_bond_harmonic * dz2;
		if(flags->eflag) vars->Uvap+=rk*dr;

		C = dx1*dx2 + dy1*dy2 + dz1*dz2;
		C /= r1*r2;
		Cs = 1/(sqrt(1.0-C*C));
		dtheta = acos(C) - 1.8242;
		tk = 55 * dtheta;
		a = -2.0 * tk * Cs;
		a11 = a*C / rsq1;
		a12 = -a / (r1*r2);
		a22 = a*C / rsq2;
		f1[0] = a11*dx1 + a12*dx2;
		f1[1] = a11*dy1 + a12*dy2;
		f1[2] = a11*dz1 + a12*dz2;
		f3[0] = a22*dx2 + a12*dx1;
		f3[1] = a22*dy2 + a12*dy1;
		f3[2] = a22*dz2 + a12*dz1;
		vapors[I].inAtoms[1].fx += f1[0];
		vapors[I].inAtoms[1].fy += f1[1];
		vapors[I].inAtoms[1].fz += f1[2];
		vapors[I].inAtoms[0].fx -= f1[0] + f3[0];
		vapors[I].inAtoms[0].fy -= f1[1] + f3[1];
		vapors[I].inAtoms[0].fz -= f1[2] + f3[2];
		vapors[I].inAtoms[2].fx += f3[0];
		vapors[I].inAtoms[2].fy += f3[1];
		vapors[I].inAtoms[2].fz += f3[2];
		if (flags->eflag) vars->Uvap += tk*dtheta;
	}
}
