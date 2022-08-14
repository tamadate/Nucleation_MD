//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


void
MD::makeDiatomicProp_in(int i){
	Molecule *gases = vars->gases.data();
    if(gases[i].inFlag==0){
        random_device seed;
        mt19937 mt(seed());
        normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
        uniform_real_distribution<double> r(0,1);

        double a=r(mt)*2*M_PI;
        double b=r(mt)*2*M_PI;
    //    double c=r(mt)*2*M_PI;

		gases[i].inAtoms=vars->atomGas(); //check name

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
					dx=gases[i].qx-ion_r[0];
					dy=gases[i].qy-ion_r[1];
					dz=gases[i].qz-ion_r[2];
					r2=dx*dx+dy*dy+dz*dz;
					double rinv=1/sqrt(r2);
					gases[i].qx-=dx*rinv;
					gases[i].qy-=dy*rinv;
					gases[i].qz-=dz*rinv;
					breakFlag=0;
				}
				if(breakFlag==0) break;
			}
		}


		int GmolSize=gases[i].inAtoms.size();
		double vy_rot = distgas(mt) *1e-5*0.1;
		double vz_rot = distgas(mt) *1e-5*0.1;
		gases[i].inAtoms[0].px=0;
		gases[i].inAtoms[0].py=vy_rot;
		gases[i].inAtoms[0].pz=0;
		gases[i].inAtoms[1].px=0;
		gases[i].inAtoms[1].py=-vy_rot;
		gases[i].inAtoms[1].pz=0;
		for (int j=0; j<GmolSize; j++){
		    //vars->ROTATION(gases[i].inAtoms[j].qx,gases[i].inAtoms[j].qy,gases[i].inAtoms[j].qz,0,a,b,gases[i].inAtoms[j].qx,gases[i].inAtoms[j].qy,gases[i].inAtoms[j].qz);
		    //vars->ROTATION(gases[i].inAtoms[j].px,gases[i].inAtoms[j].py,gases[i].inAtoms[j].pz,a,b,c,gases[i].inAtoms[j].px,gases[i].inAtoms[j].py,gases[i].inAtoms[j].pz);
			gases[i].inAtoms[j].qx+=gases[i].qx;
			gases[i].inAtoms[j].qy+=gases[i].qy;
			gases[i].inAtoms[j].qz+=gases[i].qz;
			gases[i].inAtoms[j].px+=gases[i].px;
			gases[i].inAtoms[j].py+=gases[i].py;
			gases[i].inAtoms[j].pz+=gases[i].pz;
			cout<<gases[i].inAtoms[j].px<<" "<<gases[i].inAtoms[j].py<<" "<<gases[i].inAtoms[j].pz<<" ";
		}
		cout<<endl;

    }
	vars->gases[i].inFlag=1;
}

void
MD::makeDiatomicProp_out(int i){
    Molecule *gases = vars->gases.data();
    if(gases[i].inFlag==1){
		int GmolSize=gases[i].inAtoms.size();
		for (int j=0; j<GmolSize; j++){
			gases[i].inAtoms[j].qx-=gases[i].qx;
			gases[i].inAtoms[j].qy-=gases[i].qy;
			gases[i].inAtoms[j].qz-=gases[i].qz;
			gases[i].inAtoms[j].px-=gases[i].px;
			gases[i].inAtoms[j].py-=gases[i].py;
			gases[i].inAtoms[j].pz-=gases[i].pz;
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
		int VmolSize=vapors[i].inAtoms.size();

		int gasInSize=vars->gas_in.size();
		int vaporInSize=vars->vapor_in.size();
		int breakFlag=0;
		while(breakFlag==0){
			breakFlag=1;
			for(int j=0; j<gasInSize; j++){
				double dx=gases[vars->gas_in[j]].qx-vapors[i].qx;
				double dy=gases[vars->gas_in[j]].qy-vapors[i].qy;
				double dz=gases[vars->gas_in[j]].qz-vapors[i].qz;
				double r2=dx*dx+dy*dy+dz*dz;
				if(r2<25){
					dx=vapors[i].qx-ion_r[0];
					dy=vapors[i].qy-ion_r[1];
					dz=vapors[i].qz-ion_r[2];
					r2=dx*dx+dy*dy+dz*dz;
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
					dx=vapors[i].qx-ion_r[0];
					dy=vapors[i].qy-ion_r[1];
					dz=vapors[i].qz-ion_r[2];
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


		for (int j=0; j<VmolSize; j++){
		    vars->ROTATION(vapors[i].inAtoms[j].qx,vapors[i].inAtoms[j].qy,vapors[i].inAtoms[j].qz,a,b,c,vapors[i].inAtoms[j].qx,vapors[i].inAtoms[j].qy,vapors[i].inAtoms[j].qz);
		    //vars->ROTATION(vapors[i].inAtoms[j].px,vapors[i].inAtoms[j].py,vapors[i].inAtoms[j].pz,a,b,c,vapors[i].inAtoms[j].px,vapors[i].inAtoms[j].py,vapors[i].inAtoms[j].pz);
			vapors[i].inAtoms[j].qx+=vapors[i].qx;
			vapors[i].inAtoms[j].qy+=vapors[i].qy;
			vapors[i].inAtoms[j].qz+=vapors[i].qz;
			vapors[i].inAtoms[j].px+=vapors[i].px;
			vapors[i].inAtoms[j].py+=vapors[i].py;
			vapors[i].inAtoms[j].pz+=vapors[i].pz;
		}
		if(collisionFlagVapor[i]==0){
			Ovin(i);
			collisionFlagVapor[i]=itime;
		}

    }
	vars->vapors[i].inFlag=1;
}


void
MD::makePolyatomicProp_out(int i){
    Molecule *vapors = vars->vapors.data();
    if(vapors[i].inFlag==1){
		int VmolSize=vapors[i].inAtoms.size();
		for (int j=0; j<VmolSize; j++){
			vapors[i].inAtoms[j].qx-=vapors[i].qx;
			vapors[i].inAtoms[j].qy-=vapors[i].qy;
			vapors[i].inAtoms[j].qz-=vapors[i].qz;
			vapors[i].inAtoms[j].px-=vapors[i].px;
			vapors[i].inAtoms[j].py-=vapors[i].py;
			vapors[i].inAtoms[j].pz-=vapors[i].pz;
		}
		if(collisionFlagVapor[i]!=0){
			output_vapor_collision(collisionFlagVapor[i]);
			Ovout(i);
			collisionFlagVapor[i]=0;
		}
    }
	vars->vapors[i].inFlag=0;
}
