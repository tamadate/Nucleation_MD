#include "potentialAMBER.hpp"

void
PotentialAMBER::computeBond(Variables *vars, FLAG *flags) {
	Atom *ions = vars->Molecules[0].inAtoms.data();
	Bond_type *btypes = vars->btypes.data();
	Bond *bonds=vars->Molecules[0].bonds.data();
	int bsize=vars->Molecules[0].bonds.size();
	for (int ib=0;ib<bsize;ib++) {
		int i=bonds[ib].atom1;
		int j=bonds[ib].atom2;
		int type=bonds[ib].type;
		double dx = ions[i].qx - ions[j].qx;
		double dy = ions[i].qy - ions[j].qy;
		double dz = ions[i].qz - ions[j].qz;
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r = sqrt(rsq);
		double dr = (r-btypes[type].coeff[1]);
		double rk = btypes[type].coeff[0] * dr;
		double force_bond_harmonic;
		force_bond_harmonic = -2.0*rk/r;
		ions[i].fx += force_bond_harmonic * dx;
		ions[i].fy += force_bond_harmonic * dy;
		ions[i].fz += force_bond_harmonic * dz;
		ions[j].fx -= force_bond_harmonic * dx;
		ions[j].fy -= force_bond_harmonic * dy;
		ions[j].fz -= force_bond_harmonic * dz;
		if(flags->eflag) vars->Utotal.Uion+=rk*dr;
	}
}
