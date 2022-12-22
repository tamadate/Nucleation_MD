//------------------------------------------------------------------------
#include "variables.hpp"
//------------------------------------------------------------------------
std::vector<Atom>
Variables::atomGas(int gastype) {
 	std::vector<Atom> atms;
	Atom a;
	a.type=gastype;
	a.px=a.py=a.pz=a.fx=a.fy=a.fz=a.charge=a.ix=a.iy=a.iz=0;


	if(gastype==1){
		a.id = 0;
		a.qx = 0;
		a.qy = 0;
		a.qz = 0;
		a.mass = 4.027;
		atms.push_back(a);
	}
	if(gastype==2){
		a.id = 0;
		a.qx = -0.549;
		a.qy = 0;
		a.qz = 0;
		a.mass = 14.01;
		atms.push_back(a);

		a.id = 1;
		a.qx = 0.549;
		a.qy = 0;
		a.qz = 0;
		a.mass = 14.01;
		atms.push_back(a);
	}
	if(gastype==3){
		a.id = 0;
		a.qx = 0;
		a.qy = 0;
		a.qz = 0;
		a.mass = 28.02;
		atms.push_back(a);
	}
	if(gastype==4){
		a.id = 0;
		a.qx = 0;
		a.qy = 0;
		a.qz = 0;
		a.mass = 28.02;
		atms.push_back(a);
	}
	return atms;
}
