#include "potentialAMBER.hpp"

void
PotentialAMBER::computeAngle(Variables *vars, FLAG *flags) {
	Atom *ions = vars->Molecules[0].inAtoms.data();
	Angle_type *ctypes = vars->ctypes.data();
	Angle *angles=vars->Molecules[0].angles.data();
	int asize=vars->Molecules[0].angles.size();
	for (int ian=0;ian<asize;ian++) {
		double dx1, dy1, dz1, dx2, dy2, dz2, rsq1, rsq2, r1, r2, C, Cs, dtheta, tk, a, a11, a12, a22, f1[3], f3[3];
		int i=angles[ian].atom1;
		int j=angles[ian].atom2;
		int k=angles[ian].atom3;
		int type=angles[ian].type;
		dx1 = ions[i].qx - ions[j].qx;
		dy1 = ions[i].qy - ions[j].qy;
		dz1 = ions[i].qz - ions[j].qz;
		rsq1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
		r1 = sqrt(rsq1);
		dx2 = ions[k].qx - ions[j].qx;
		dy2 = ions[k].qy - ions[j].qy;
		dz2 = ions[k].qz - ions[j].qz;
		rsq2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
		r2 = sqrt(rsq2);
		C = dx1*dx2 + dy1*dy2 + dz1*dz2;
		C /= r1*r2;
		Cs = 1/(sqrt(1.0-C*C));
		dtheta = acos(C) - ctypes[type].coeff[1];
		tk = ctypes[type].coeff[0] * dtheta;
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
  	ions[i].fx += f1[0];
		ions[i].fy += f1[1];
		ions[i].fz += f1[2];
  	ions[j].fx -= f1[0] + f3[0];
		ions[j].fy -= f1[1] + f3[1];
		ions[j].fz -= f1[2] + f3[2];
		ions[k].fx += f3[0];
		ions[k].fy += f3[1];
		ions[k].fz += f3[2];
    if (flags->eflag) vars->Utotal.Uion+= tk*dtheta;
	}
}
