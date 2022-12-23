//------------------------------------------------------------------------
#include "../md.hpp"
//------------------------------------------------------------------------

void
MD::makeDiatomicProp_in(Molecule &gasOut){
	random_device seed;
  mt19937 mt(seed());
  normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	uniform_real_distribution<double> r(0,1);
	std::vector<Atom> at=vars->atomGas(gastype);
	gasOut.mass=pp->Mgas;
	double a=r(mt)*2*M_PI;
	double b=r(mt)*2*M_PI;
	double c=r(mt)*2*M_PI;
	double vy_rot = distgas(mt) *1e-5;
	double vz_rot = distgas(mt) *1e-5;
	if(gastype==2){
		at[0].px=0;
		at[0].py=vy_rot;
		at[0].pz=vz_rot;
		at[1].px=0;
		at[1].py=-vy_rot;
		at[1].pz=-vz_rot;
	}
	for (auto &ag : at){
    vars->ROTATION(ag.qx,ag.qy,ag.qz,a,b,c,ag.qx,ag.qy,ag.qz);
		vars->ROTATION(ag.px,ag.py,ag.pz,a,b,c,ag.px,ag.py,ag.pz);
		ag.qx+=gasOut.qx;
		ag.qy+=gasOut.qy;
		ag.qz+=gasOut.qz;
		ag.px+=gasOut.px;
		ag.py+=gasOut.py;
		ag.pz+=gasOut.pz;
  }
	gasOut.inAtoms=at;
}


void
MD::makePolyatomicProp_in(Molecule &vapOut){
  random_device seed;
  mt19937 mt(seed());
  default_random_engine engine(seed());
  normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mvapor));
  uniform_real_distribution<double> r(0,1);

  double a=r(mt)*2*M_PI;
  double b=r(mt)*2*M_PI;
  double c=r(mt)*2*M_PI;
	std::vector<Atom> at=vars->atomVapor;
	vapOut.mass=pp->mvapor;

	for (auto &av : at){
    vars->ROTATION(av.qx,av.qy,av.qz,a,b,c,av.qx,av.qy,av.qz);
    vars->ROTATION(av.px,av.py,av.pz,a,b,c,av.px,av.py,av.pz);
		av.qx+=vapOut.qx;
		av.qy+=vapOut.qy;
		av.qz+=vapOut.qz;
		av.px+=vapOut.px;
		av.py+=vapOut.py;
		av.pz+=vapOut.pz;
	}
	totalVaporIn++;
	vapOut.id=totalVaporIn;
	vapOut.inAtoms=at;

	sprintf(filepath, "vapor_in_%d.dat", int(calculation_number));
	FILE*f=fopen(filepath, "a");
	Molecule *ion=vars->Molecules.data();
	fprintf(f, "%d %e %e %e %e %e %e %e\n", totalVaporIn,itime*dt,vapOut.qx-ion[0].qx,vapOut.qy-ion[0].qy,vapOut.qz-ion[0].qz,\
	vapOut.px-ion[0].px,vapOut.py-ion[0].py,vapOut.pz-ion[0].pz);
	fclose(f);
}
