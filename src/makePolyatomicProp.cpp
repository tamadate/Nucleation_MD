//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


void 
MD::makeDiatomicProp_in(int i){
    
    if(vars->gases[i+Nof_around_gas].mass==0){
        random_device seed;
        mt19937 mt(seed());
        default_random_engine engine(seed());
        normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
        uniform_real_distribution<double> r(0,1);
        
        const double NNbond=1.098;
        const double NNbond_half=NNbond*0.5;
        //const double I1=0.5*pp->mgas*NNbond_half*NNbond_half;
        double a=r(mt)*2*M_PI;
        double b=r(mt)*2*M_PI;
        //double cosb=1-2*r(mt);
        //double b=acos(cosb);
        //double c=r(mt)*2*M_PI;
        //double wsq=-2*kb*T/I1/log(1-r(mt));
        //double w2=sqrt(wsq/(1+tan(c)*tan(c)))*1e-5;
        //double w1=tan(c)*w2;
        //cout<<sqrt(wsq)<<" "<<w1<<" "<<w2<<" "<<endl;
        
        double x1=-NNbond_half;
        double y1=0;
        double z1=0;
        double x2=NNbond_half;
        double y2=0;
        double z2=0;
        
        double vx1=0;
        double vy1=distgas(engine)*1e-5;
        double vz1=distgas(engine)*1e-5;
        double vx2=0;
        double vy2=-vy1;
        double vz2=-vz1;
        
        vars->ROTATION(x1,y1,z1,0,a,b,x1,y1,z1);
        vars->ROTATION(x2,y2,z2,0,a,b,x2,y2,z2);
        vars->ROTATION(vx1,vy1,vz1,0,a,b,vx1,vy1,vz1);
        vars->ROTATION(vx2,vy2,vz2,0,a,b,vx2,vy2,vz2);

        
        //  set positions
        //  index i have a center of gas molecule and i+Nof... have no information before this operation
        vars->gases[i+Nof_around_gas].qx=vars->gases[i].qx+x2;
        vars->gases[i+Nof_around_gas].qy=vars->gases[i].qy+y2;
        vars->gases[i+Nof_around_gas].qz=vars->gases[i].qz+z2;
        vars->gases[i].qx+=x1;
        vars->gases[i].qy+=y1;
        vars->gases[i].qz+=z1;
        //  set velocities
        vars->gases[i+Nof_around_gas].px=vars->gases[i].px+vx2;
        vars->gases[i+Nof_around_gas].py=vars->gases[i].py+vy2;
        vars->gases[i+Nof_around_gas].pz=vars->gases[i].pz+vz2;
        vars->gases[i].px+=vx1;
        vars->gases[i].py+=vy1;
        vars->gases[i].pz+=vz1;
        
        vars->gases[i].mass=14.01;
        vars->gases[i+Nof_around_gas].mass=14.01;
    }
}

void 
MD::makeDiatomicProp_out(int i){
    if(vars->gases[i+Nof_around_gas].mass!=0){
        //  averaged position
        vars->gases[i].qx=(vars->gases[i].qx+vars->gases[i+Nof_around_gas].qx)*0.5;
        vars->gases[i].qy=(vars->gases[i].qy+vars->gases[i+Nof_around_gas].qy)*0.5;
        vars->gases[i].qz=(vars->gases[i].qz+vars->gases[i+Nof_around_gas].qz)*0.5;
        //  averaged velocity
        vars->gases[i].px=(vars->gases[i].px+vars->gases[i+Nof_around_gas].px)*0.5;
        vars->gases[i].py=(vars->gases[i].py+vars->gases[i+Nof_around_gas].py)*0.5;
        vars->gases[i].pz=(vars->gases[i].pz+vars->gases[i+Nof_around_gas].pz)*0.5;
        
        vars->gases[i].mass=28.02;
        vars->gases[i+Nof_around_gas].mass=0;
        
        vars->gases[i+Nof_around_gas].px=0;
        vars->gases[i+Nof_around_gas].py=0;
        vars->gases[i+Nof_around_gas].pz=0;

		if(collisionFlagGas[i]!=0){
			output_gas_collision(collisionFlagGas[i]);
			collisionFlagGas[i]=0;
		}
    }
}

void 
MD::makePolyatomicProp_in(int i){
    Molecule *vapors = vars->vapors.data();
	Atom *gases = vars->gases.data();
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

