//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------

void
MD::makeDiatomicProp_in(Molecule_out &gasOut){
	random_device seed;
  mt19937 mt(seed());
  normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	uniform_real_distribution<double> r(0,1);
	int breakFlag=0;
	while(breakFlag==0){
		breakFlag=1;
		for(auto &gin : vars->effectiveIn[1]){
			double dx=gin.qx-gasOut.qx;
			double dy=gin.qy-gasOut.qy;
			double dz=gin.qz-gasOut.qz;
			double r2=dx*dx+dy*dy+dz*dz;
			if(r2<25){
				dx=gasOut.qx-vars->effectiveIn[0][0].qx;
				dy=gasOut.qy-vars->effectiveIn[0][0].qy;
				dz=gasOut.qz-vars->effectiveIn[0][0].qz;
				r2=dx*dx+dy*dy+dz*dz;
				double rinv=1/sqrt(r2);
				gasOut.qx-=dx*rinv;
				gasOut.qy-=dy*rinv;
				gasOut.qz-=dz*rinv;
				breakFlag=0;
			}
			if(breakFlag==0) break;
		}
	}
	Molecule mol;
	mol.inAtoms=vars->atomGas(gastype);
	mol.mass=pp->Mgas;
	for (auto &b : mol.inAtoms){
		for (int thread=0;thread<Nth;thread++){
			b.fxMP.push_back(0);
			b.fyMP.push_back(0);
			b.fzMP.push_back(0);
		}
	}
	double a=r(mt)*2*M_PI;
	double b=r(mt)*2*M_PI;
	double c=r(mt)*2*M_PI;
	double vy_rot = distgas(mt) *1e-5;
	double vz_rot = distgas(mt) *1e-5;
	if(gastype==2){
		mol.inAtoms[0].px=0;
		mol.inAtoms[0].py=vy_rot;
		mol.inAtoms[0].pz=vz_rot;
		mol.inAtoms[1].px=0;
		mol.inAtoms[1].py=-vy_rot;
		mol.inAtoms[1].pz=-vz_rot;
	}
	for (auto &ag : mol.inAtoms){
    vars->ROTATION(ag.qx,ag.qy,ag.qz,a,b,c,ag.qx,ag.qy,ag.qz);
		vars->ROTATION(ag.px,ag.py,ag.pz,a,b,c,ag.px,ag.py,ag.pz);
		ag.qx+=gasOut.qx;
		ag.qy+=gasOut.qy;
		ag.qz+=gasOut.qz;
		ag.px+=gasOut.px;
		ag.py+=gasOut.py;
		ag.pz+=gasOut.pz;
  }
	mol.qx=gasOut.qx;
	mol.qy=gasOut.qy;
	mol.qz=gasOut.qz;
	mol.px=gasOut.px;
	mol.py=gasOut.py;
	mol.pz=gasOut.pz;
	vars->effectiveIn[1].push_back(mol);
}


void
MD::makePolyatomicProp_in(Molecule_out &vapOut){
  random_device seed;
  mt19937 mt(seed());
  default_random_engine engine(seed());
  normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mvapor));
  uniform_real_distribution<double> r(0,1);

  double a=r(mt)*2*M_PI;
  double b=r(mt)*2*M_PI;
  double c=r(mt)*2*M_PI;

	int breakFlag=0;
	while(breakFlag==0){
		breakFlag=1;
		for(auto &gin : vars->effectiveIn[1]){
			double dx=gin.qx-vapOut.qx;
			double dy=gin.qy-vapOut.qy;
			double dz=gin.qz-vapOut.qz;
			double r2=dx*dx+dy*dy+dz*dz;
			if(r2<25){
				dx=vapOut.qx-vars->effectiveIn[0][0].qx;
				dy=vapOut.qy-vars->effectiveIn[0][0].qy;
				dz=vapOut.qz-vars->effectiveIn[0][0].qz;
				r2=dx*dx+dy*dy+dz*dz;
				double rinv=1/sqrt(r2);
				vapOut.qx-=dx*rinv;
				vapOut.qy-=dy*rinv;
				vapOut.qz-=dz*rinv;
				breakFlag=0;
			}
			if(breakFlag==0) break;
		}

		for(auto &vin : vars->effectiveIn[2]){
			double dx=vin.qx-vapOut.qx;
			double dy=vin.qy-vapOut.qy;
			double dz=vin.qz-vapOut.qz;
			double r2=dx*dx+dy*dy+dz*dz;
			if(r2<25){
				dx=vapOut.qx-vars->effectiveIn[0][0].qx;
				dy=vapOut.qy-vars->effectiveIn[0][0].qy;
				dz=vapOut.qz-vars->effectiveIn[0][0].qz;
				r2=dx*dx+dy*dy+dz*dz;
				double rinv=1/sqrt(r2);
				vapOut.qx-=dx*rinv;
				vapOut.qy-=dy*rinv;
				vapOut.qz-=dz*rinv;
				breakFlag=0;
			}
		}
	}

	Molecule mol;
	mol.inAtoms=vars->atomVapor;
	for (auto &b : mol.inAtoms){
		for (int thread=0;thread<Nth;thread++){
			b.fxMP.push_back(0);
			b.fyMP.push_back(0);
			b.fzMP.push_back(0);
		}
	}
	mol.bonds=vars->bonds_v;
	mol.angles=vars->angles_v;
	mol.dihedrals=vars->dihedrals_v;

	for (auto &av : mol.inAtoms){
    vars->ROTATION(av.qx,av.qy,av.qz,a,b,c,av.qx,av.qy,av.qz);
    vars->ROTATION(av.px,av.py,av.pz,a,b,c,av.px,av.py,av.pz);
		av.qx+=vapOut.qx;
		av.qy+=vapOut.qy;
		av.qz+=vapOut.qz;
		av.px+=vapOut.px;
		av.py+=vapOut.py;
		av.pz+=vapOut.pz;
	}
	mol.qx=vapOut.qx;
	mol.qy=vapOut.qy;
	mol.qz=vapOut.qz;
	mol.px=vapOut.px;
	mol.py=vapOut.py;
	mol.pz=vapOut.pz;
	totalVaporIn++;
	mol.id=totalVaporIn;

	sprintf(filepath, "vapor_in_%d.dat", int(calculation_number));
	FILE*f=fopen(filepath, "a");
	Molecule *ion=vars->effectiveIn[0].data();
	fprintf(f, "%d %e %e %e %e %e %e %e\n", totalVaporIn,itime*dt,vapOut.qx-ion[0].qx,vapOut.qy-ion[0].qy,vapOut.qz-ion[0].qz,\
	vapOut.px-ion[0].px,vapOut.py-ion[0].py,vapOut.pz-ion[0].pz);
	fclose(f);

	vars->effectiveIn[2].push_back(mol);
}
