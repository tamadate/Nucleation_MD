//------------------------------------------------------------------------
#include "variables.hpp"
//------------------------------------------------------------------------
void
Variables::add_atype(int type, double mass, double coeff1, double coeff2) {
 	Atom_type at;
	at.type = type;
	at.mass = mass; 
	at.coeff1 = coeff1;
	at.coeff2 = coeff2;
	atypes.push_back(at);
}
void
Variables::add_ions(int id, int type, double x, double y, double z, double vx, double vy, double vz, double fx, double fy, double fz, double charge, double mass) {
 	Ion a;
	a.id = id;
	a.type = type; 
	a.qx = x;
	a.qy = y;
	a.qz = z;
	a.px = vx;
	a.py = vy;
	a.pz = vz;
	a.mass = mass;
	a.fx = fx;
	a.fy = fy;
	a.fz = fz;
	a.charge = charge;
	a.ix=0;
	a.iy=0;
	a.iz=0;
	ions.push_back(a);
}
void
Variables::add_gases(int id, int type, double x, double y, double z, double vx, double vy, double vz, double fx, double fy, double fz, double charge, double mass) {
 	Gas a;
	a.id = id;
	a.type = type; 
	a.qx = x;
	a.qy = y;
	a.qz = z;
	a.px = vx;
	a.py = vy;
	a.pz = vz;
	a.mass = mass;
	a.fx = fx;
	a.fy = fy;
	a.fz = fz;
	a.charge = charge;
	a.ix=0;
	a.iy=0;
	a.iz=0;
	gases.push_back(a);
}
void
Variables::add_vapors(int id, int type, double x, double y, double z, double vx, double vy, double vz, double mass) {
 	Molecule a;
	a.qx = x;
	a.qy = y;
	a.qz = z;
	a.px = vx;
	a.py = vy;
	a.pz = vz;
	a.ix=0;
	a.iy=0;
	a.iz=0;
	Atom b;
	a.inAtoms=atomVapor;
	a.bonds=bondMeOH();
	a.angles=angleMeOH();
	a.dihedrals=dihedralMeOH();
	a.inFlag=0;
	vapors.push_back(a);
}
void
Variables::add_bonds(int atom1, int atom2, int type) {
 	Bond b;
	b.atom1 = atom1;
	b.atom2 = atom2; 
	b.type = type;
	bonds.push_back(b);
}
void
Variables::add_angles(int atom1, int atom2, int atom3, int type) {
 	Angle c;
	c.atom1 = atom1;
	c.atom2 = atom2; 
	c.atom3 = atom3;
	c.type = type;
	angles.push_back(c);
}
void
Variables::add_dihedrals(int atom1, int atom2, int atom3, int atom4, int type) {
 	Dihedral d;
	d.atom1 = atom1;
	d.atom2 = atom2;
	d.atom3 = atom3;
	d.atom4 = atom4;
	d.type = type;
	dihedrals.push_back(d);
}
