//------------------------------------------------------------------------
#include "../variables.hpp"

double
Variables::distFromIon(Molecule mol) {
	double dx=mol.qx-Molecules[0].qx;
  double dy=mol.qy-Molecules[0].qy;
  double dz=mol.qz-Molecules[0].qz;
  return dx*dx+dy*dy+dz*dz;
}
