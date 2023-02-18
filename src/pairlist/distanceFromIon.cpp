#include "../md.hpp"


double
MD::distFromIonCenter(Molecule &mol, double &dx, double &dy, double &dz) {
	dx = mol.qx - vars->ion_r[0];
	dy = mol.qy - vars->ion_r[1];
	dz = mol.qz - vars->ion_r[2];
	return dx*dx+dy*dy+dz*dz;
}
