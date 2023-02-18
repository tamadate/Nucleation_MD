#include "../md.hpp"


void
MD::makeDiatomicProp_in(int i){
	Molecule *gases = vars->gases.data();
  if(gases[i].inFlag==0){
  	random_device seed;
    mt19937 mt(seed());
    normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
		uniform_real_distribution<double> r(0,1);

		gases[i].inAtoms=vars->atomGas(gastype); //check name

		int gasInSize=vars->gas_in.size();
		int breakFlag=0;
		while(breakFlag==0){
			breakFlag=1;
			for(int j=0; j<gasInSize-1; j++){
				double dx=gases[vars->gas_in[j]].qx-gases[i].qx;
				double dy=gases[vars->gas_in[j]].qy-gases[i].qy;
				double dz=gases[vars->gas_in[j]].qz-gases[i].qz;
				double r2=dx*dx+dy*dy+dz*dz;
				if(r2<25){
					r2=distFromIonCenter(gases[i],dx,dy,dz);
					double rinv=1/sqrt(r2);
					gases[i].qx-=dx*rinv;
					gases[i].qy-=dy*rinv;
					gases[i].qz-=dz*rinv;
					breakFlag=0;
				}
				if(breakFlag==0) break;
			}
		}

		double a=r(mt)*2*M_PI;
		double b=r(mt)*2*M_PI;
		double c=r(mt)*2*M_PI;
		double vy_rot = distgas(mt) *1e-5;
		double vz_rot = distgas(mt) *1e-5;
		if(gastype==2){
			gases[i].inAtoms[0].px=0;
			gases[i].inAtoms[0].py=vy_rot;
			gases[i].inAtoms[0].pz=vz_rot;
			gases[i].inAtoms[1].px=0;
			gases[i].inAtoms[1].py=-vy_rot;
			gases[i].inAtoms[1].pz=-vz_rot;
		}
		for (auto &ag : gases[i].inAtoms){
	    vars->ROTATION(ag.qx,ag.qy,ag.qz,a,b,c,ag.qx,ag.qy,ag.qz);
			vars->ROTATION(ag.px,ag.py,ag.pz,a,b,c,ag.px,ag.py,ag.pz);
			ag.qx+=gases[i].qx;
			ag.qy+=gases[i].qy;
			ag.qz+=gases[i].qz;
			ag.px+=gases[i].px;
			ag.py+=gases[i].py;
			ag.pz+=gases[i].pz;
		}
  }
	vars->gases[i].inFlag=1;
}

void
MD::makeDiatomicProp_out(int i){
  Molecule *gases = vars->gases.data();
  if(gases[i].inFlag==1){
		for (auto &a : gases[i].inAtoms){
			a.qx-=gases[i].qx;
			a.qy-=gases[i].qy;
			a.qz-=gases[i].qz;
			a.px-=gases[i].px;
			a.py-=gases[i].py;
			a.pz-=gases[i].pz;
		}
  }
	vars->gases[i].inFlag=0;
}

void
MD::makePolyatomicProp_in(int i){
  	Molecule *vapors = vars->vapors.data();
	Molecule *gases = vars->gases.data();
  	if(vapors[i].inFlag==0){
		random_device seed;
		mt19937 mt(seed());
		default_random_engine engine(seed());
		normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mvapor));
		uniform_real_distribution<double> r(0,1);

		double a=r(mt)*2*M_PI;
		double b=r(mt)*2*M_PI;
		double c=r(mt)*2*M_PI;

		vapors[i].inAtoms=vars->atomVapor; //check name

		int vaporInSize=vars->vapor_in.size();
		int breakFlag=0;
		while(breakFlag==0){
			breakFlag=1;
			for(auto j : vars->gas_in){
				double dx=gases[j].qx-vapors[i].qx;
				double dy=gases[j].qy-vapors[i].qy;
				double dz=gases[j].qz-vapors[i].qz;
				double r2=dx*dx+dy*dy+dz*dz;
				if(r2<25){
					r2=distFromIonCenter(vapors[i],dx,dy,dz);
					double rinv=1/sqrt(r2);
					vapors[i].qx-=dx*rinv;
					vapors[i].qy-=dy*rinv;
					vapors[i].qz-=dz*rinv;
					breakFlag=0;
				}
				if(breakFlag==0) break;
			}


			for(int j=0; j<vaporInSize-1; j++){
				double dx=vapors[vars->vapor_in[j]].qx-vapors[i].qx;
				double dy=vapors[vars->vapor_in[j]].qy-vapors[i].qy;
				double dz=vapors[vars->vapor_in[j]].qz-vapors[i].qz;
				double r2=dx*dx+dy*dy+dz*dz;
				if(r2<25){
					dx=vapors[i].qx-vars->ion_r[0];
					dy=vapors[i].qy-vars->ion_r[1];
					dz=vapors[i].qz-vars->ion_r[2];
					r2=dx*dx+dy*dy+dz*dz;
					double rinv=1/sqrt(r2);
					vapors[i].qx-=dx*rinv;
					vapors[i].qy-=dy*rinv;
					vapors[i].qz-=dz*rinv;
					breakFlag=0;
				}
				if(breakFlag==0) break;
			}
		}


		for (auto &av : vapors[i].inAtoms){
		vars->ROTATION(av.qx,av.qy,av.qz,a,b,c,av.qx,av.qy,av.qz);
		//vars->ROTATION(av.px,av.py,av.pz,a,b,c,av.px,av.py,av.pz);
			av.qx+=vapors[i].qx;
			av.qy+=vapors[i].qy;
			av.qz+=vapors[i].qz;
			av.px+=vapors[i].px;
			av.py+=vapors[i].py;
			av.pz+=vapors[i].pz;
		}

  	}
	vars->vapors[i].inFlag=1;
}


void
MD::makePolyatomicProp_out(int i){
  	Molecule *vapors = vars->vapors.data();
  	if(vapors[i].inFlag==1){
		for (auto &av : vapors[i].inAtoms){
			av.qx-=vapors[i].qx;
			av.qy-=vapors[i].qy;
			av.qz-=vapors[i].qz;
			av.px-=vapors[i].px;
			av.py-=vapors[i].py;
			av.pz-=vapors[i].pz;
		}
		if(collisionFlagVapor[i]!=0){
			obs->Ovout(i,itime*dt);
			collisionFlagVapor[i]=0;
		}
  	}
	vars->vapors[i].inFlag=0;
}
