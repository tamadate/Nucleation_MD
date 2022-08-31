//------------------------------------------------------------------------
#include "variables.hpp"
//------------------------------------------------------------------------
std::vector<Atom>
Variables::makeAtomTIP3P(void) {
 	std::vector<Atom> atms;
	Atom a;
	a.id = 0;
	a.type = 9;
	a.qx = 0;
	a.qy = 0;
	a.qz = 0;
	a.px = 0;
	a.py = 0;
	a.pz = 0;
	a.mass = 15.999;
	a.fx = 0;
	a.fy = 0;
	a.fz = 0;
	a.charge = -0.834;
	a.ix=0;
	a.iy=0;
	a.iz=0;
  
  for (int thread=0;thread<Nth;thread++){
    a.fxMP.push_back(0);
    a.fyMP.push_back(0);
    a.fzMP.push_back(0);
  }
	atms.push_back(a);

	a.id = 1;
	a.type = 10;
	a.qx = 0.75695;
	a.qy = -0.58588;
	a.qz = 0;
	a.mass = 1.008;
	a.charge = 0.417;
	atms.push_back(a);

	a.id = 2;
	a.type = 10;
	a.qx = -0.75695;
	a.qy = -0.58588;
	a.qz = 0;
	a.mass = 1.008;
	a.charge = 0.417;
	atms.push_back(a);

	return atms;
}
