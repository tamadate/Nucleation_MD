#include "../md.hpp"


double
MD::distFromIonCenter(Molecule &mol, double &dx, double &dy, double &dz) {
	dx = mol.qx - vars->IonX[0];
	dy = mol.qy - vars->IonX[1];
	dz = mol.qz - vars->IonX[2];
	return dx*dx+dy*dy+dz*dz;
}


double
MD::distFromIonPreCenter(Molecule &mol, double &dx, double &dy, double &dz) {
	dx = mol.qx - vars->preIonX[0];
	dy = mol.qy - vars->preIonX[1];
	dz = mol.qz - vars->preIonX[2];
	return dx*dx+dy*dy+dz*dz;
}
