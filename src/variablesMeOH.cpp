//------------------------------------------------------------------------
#include "variables.hpp"
//------------------------------------------------------------------------
std::vector<Atom>
Variables::makeAtomMeOH(void) {
 	std::vector<Atom> atms;
	Atom a;
	a.id = 0;
	a.type = 5; 
	a.qx = -0.71259375;
	a.qy = 0.06896875;
	a.qz = 0.03515625;
	a.px = 0;
	a.py = 0;
	a.pz = 0;
	a.mass = 12;
	a.fx = 0;
	a.fy = 0;
	a.fz = 0;
	a.charge = 0.145;
	a.ix=0;
	a.iy=0;
	a.iz=0;
	atms.push_back(a);

	a.id = 1;
	a.type = 6; 
	a.qx = 0.67140625;
	a.qy = -0.10003125;
	a.qz = -0.07584375;
	a.mass = 16;
	a.charge = -0.683;
	atms.push_back(a);

	a.id = 2;
	a.type = 8; 
	a.qx = 1.07240625;
	a.qy = 0.46096875;
	a.qz = 0.63715625;
	a.mass = 1;
	a.charge = 0.418;
	atms.push_back(a);

	a.id = 3;
	a.type = 7; 
	a.qx = -1.21259375;
	a.qy = -0.54203125;
	a.qz = -0.74484375;
	a.mass = 1;
	a.charge = 0.04;
	atms.push_back(a);

	a.id = 4;
	a.type = 7; 
	a.qx = -0.98559375;
	a.qy = 1.13396875;
	a.qz = -0.12984375;
	a.mass = 1;
	a.charge = 0.04;
	atms.push_back(a);

	a.id = 5;
	a.type = 7; 
	a.qx = -1.06559375;
	a.qy = -0.28003125;
	a.qz = 1.02915625;
	a.mass = 1;
	a.charge = 0.04;
	atms.push_back(a);

	return atms;
}

std::vector<Bond>
Variables::bondMeOH(void) {
 	std::vector<Bond> bondsMeOH;
	Bond b;
	b.atom1 = 0;
	b.atom2 = 1;
	b.type = 1;
	bondsMeOH.push_back(b);
	b.atom1 = 0;
	b.atom2 = 3;
	b.type = 0;
	bondsMeOH.push_back(b);
	b.atom1 = 0;
	b.atom2 = 4;
	b.type = 0;
	bondsMeOH.push_back(b);
	b.atom1 = 0;
	b.atom2 = 5;
	b.type = 0;
	bondsMeOH.push_back(b);
	b.atom1 = 1;
	b.atom2 = 2;
	b.type = 2;
	bondsMeOH.push_back(b);

	return bondsMeOH;
}

std::vector<Angle>
Variables::angleMeOH(void) {
 	std::vector<Angle> anglesMeOH;
	Angle c;
	c.atom1 = 1;
	c.atom2 = 0;
	c.atom3 = 3;
	c.type = 1;
	anglesMeOH.push_back(c);
	c.atom1 = 1;
	c.atom2 = 0;
	c.atom3 = 4;
	c.type = 1;
	anglesMeOH.push_back(c);
	c.atom1 = 1;
	c.atom2 = 0;
	c.atom3 = 5;
	c.type = 1;
	anglesMeOH.push_back(c);
	c.atom1 = 3;
	c.atom2 = 0;
	c.atom3 = 4;
	c.type = 0;
	anglesMeOH.push_back(c);
	c.atom1 = 3;
	c.atom2 = 0;
	c.atom3 = 5;
	c.type = 0;
	anglesMeOH.push_back(c);
	c.atom1 = 4;
	c.atom2 = 0;
	c.atom3 = 5;
	c.type = 0;
	anglesMeOH.push_back(c);
	c.atom1 = 0;
	c.atom2 = 1;
	c.atom3 = 2;
	c.type = 2;
	anglesMeOH.push_back(c);

	return anglesMeOH;
}

std::vector<Dihedral>
Variables::dihedralMeOH(void) {
 	std::vector<Dihedral> dihedralsMeOH;
	Dihedral d;
	d.atom1 = 2;
	d.atom2 = 1;
	d.atom3 = 0;
	d.atom4 = 3;
	d.type = 1;
	dihedralsMeOH.push_back(d);
	d.atom1 = 2;
	d.atom2 = 1;
	d.atom3 = 0;
	d.atom4 = 4;
	d.type = 1;
	dihedralsMeOH.push_back(d);
	d.atom1 = 2;
	d.atom2 = 1;
	d.atom3 = 0;
	d.atom4 = 5;
	d.type = 1;
	dihedralsMeOH.push_back(d);

	return dihedralsMeOH;
}


