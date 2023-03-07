#include "md.hpp"
#include <random>
#include <algorithm>
#include <dirent.h>
#include <vector>
/*########################################################################################

-----Periodic conditions-----

#######################################################################################*/

/**************************periodic**********************************/
void
adjust_periodic(double &dx, double &dy, double &dz, double d_size) {
	const double LH = d_size * 0.5;
	if (dx < -LH)dx += d_size;
	if (dx > LH) dx -= d_size;
	if (dy < -LH)dy += d_size;
	if (dy > LH) dy -= d_size;
	if (dz < -LH)dz += d_size;
	if (dz > LH) dz -= d_size;
}
/************************periodec condition for in gas***************************/
void
MD::periodic(void) {
	Molecule *gases = vars->gases.data();
	int gis=vars->gas_in.size();
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;
	for (int k=0;k<gis;k++) {
		int i=vars->gas_in[k];
		flag=flagx=flagy=flagz=0;
		if (gases[i].qx < ion_r[0]-HL) gases[i].qx += d_size, flagx--, flag++;
		if (gases[i].qy < ion_r[1]-HL) gases[i].qy += d_size, flagy--, flag++;
		if (gases[i].qz < ion_r[2]-HL) gases[i].qz += d_size, flagz--, flag++;
		if (gases[i].qx > ion_r[0]+HL) gases[i].qx -= d_size, flagx++, flag++;
		if (gases[i].qy > ion_r[1]+HL) gases[i].qy -= d_size, flagy++, flag++;
		if (gases[i].qz > ion_r[2]+HL) gases[i].qz -= d_size, flagz++, flag++;
	}
}

/************************periodec condition for in gas***************************/


void
MD::boundary_scaling_gas_move(void){

	Molecule *gases = vars->gases.data();
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;

	for (auto &i : vars->gas_out){
		flag=flagx=flagy=flagz=0;
		if (gases[i].qx < pre_ion[0]-HL) gases[i].qx += d_size, flagx--, flag++;
		if (gases[i].qy < pre_ion[1]-HL) gases[i].qy += d_size, flagy--, flag++;
		if (gases[i].qz < pre_ion[2]-HL) gases[i].qz += d_size, flagz--, flag++;
		if (gases[i].qx > pre_ion[0]+HL) gases[i].qx -= d_size, flagx++, flag++;
		if (gases[i].qy > pre_ion[1]+HL) gases[i].qy -= d_size, flagy++, flag++;
		if (gases[i].qz > pre_ion[2]+HL) gases[i].qz -= d_size, flagz++, flag++;
		if (flag>0) {
			if(mbdist->number>mbdist->vflux.size()*0.9) {mbdist->makeWeightedMB(pp->cgas,pp->mgas,T);}
		    vx=gases[i].px*gases[i].mass;
		    vy=gases[i].py*gases[i].mass;
		    vz=gases[i].pz*gases[i].mass;
				double mgasInv=1/pp->m_gas;
				vx = vx*mgasInv;
				vy = vy*mgasInv;
				vz = vz*mgasInv;
				v2 = vx*vx+vy*vy+vz*vz;
				v = sqrt(v2);
				vMB = mbdist->vflux[mbdist->number]*1e-5;
				mod_factor= vMB/v;
		    vx = gases[i].px;
		    vy = gases[i].py;
		    vz = gases[i].pz;
		    v2 = vx*vx+vy*vy+vz*vz;
		    v = sqrt(v2);
		    vMB = mbdist->vflux[mbdist->number]*1e-5;
		    mod_factor= vMB/v;
		    gases[i].px = vx * mod_factor;
		    gases[i].py = vy * mod_factor;
		    gases[i].pz = vz * mod_factor;
		    mbdist->number++;
		}
	}
}




void
MD::boundary_scaling_vapor_move(void){
	Molecule *vapors = vars->vapors.data();
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;

	for (auto &i : vars->vapor_out){
		flag=flagx=flagy=flagz=0;
		if (vapors[i].qx < pre_ion[0]-HL) vapors[i].qx += d_size, flagx--, flag++;
		if (vapors[i].qy < pre_ion[1]-HL) vapors[i].qy += d_size, flagy--, flag++;
		if (vapors[i].qz < pre_ion[2]-HL) vapors[i].qz += d_size, flagz--, flag++;
		if (vapors[i].qx > pre_ion[0]+HL) vapors[i].qx -= d_size, flagx++, flag++;
		if (vapors[i].qy > pre_ion[1]+HL) vapors[i].qy -= d_size, flagy++, flag++;
		if (vapors[i].qz > pre_ion[2]+HL) vapors[i].qz -= d_size, flagz++, flag++;
		if (flag>0) {
			if(mbdistV->number>mbdistV->vflux.size()*0.9) {mbdistV->makeWeightedMB(pp->cvapor,pp->mvapor,T);}
		    vx=vapors[i].px;
		    vy=vapors[i].py;
		    vz=vapors[i].pz;
			v2 = vx*vx+vy*vy+vz*vz;
			v = sqrt(v2);
			vMB = mbdistV->vflux[mbdistV->number]*1e-5;
			mod_factor= vMB/v;
		    vapors[i].px = vx * mod_factor;
		    vapors[i].py = vy * mod_factor;
		    vapors[i].pz = vz * mod_factor;
		    mbdistV->number++;
		}
	}
}

void
MD::boundary_scaling_ion_move(void){
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;
	random_device seed;
	mt19937 mt(seed());

	Molecule *gases = vars->gases.data();
	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));

	for (auto &i : vars->gas_out){
		flag=flagx=flagy=flagz=0;
		if (gases[i].qx < ion_r[0]-HL) gases[i].qx += d_size, flagx--, flag++;
		if (gases[i].qy < ion_r[1]-HL) gases[i].qy += d_size, flagy--, flag++;
		if (gases[i].qz < ion_r[2]-HL) gases[i].qz += d_size, flagz--, flag++;
		if (gases[i].qx > ion_r[0]+HL) gases[i].qx -= d_size, flagx++, flag++;
		if (gases[i].qy > ion_r[1]+HL) gases[i].qy -= d_size, flagy++, flag++;
		if (gases[i].qz > ion_r[2]+HL) gases[i].qz -= d_size, flagz++, flag++;
		if (flag>0) {
			gases[i].px = distgas(mt) *1e-5;
			gases[i].py = distgas(mt) *1e-5;
			gases[i].pz = distgas(mt) *1e-5;
		}
	}

	Molecule *vapors = vars->vapors.data();
	normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));

	for (auto &i : vars->vapor_out){
		flag=flagx=flagy=flagz=0;
		if (vapors[i].qx < ion_r[0]-HL) vapors[i].qx += d_size, flag++;
		if (vapors[i].qy < ion_r[1]-HL) vapors[i].qy += d_size, flag++;
		if (vapors[i].qz < ion_r[2]-HL) vapors[i].qz += d_size, flag++;
		if (vapors[i].qx > ion_r[0]+HL) vapors[i].qx -= d_size, flag++;
		if (vapors[i].qy > ion_r[1]+HL) vapors[i].qy -= d_size, flag++;
		if (vapors[i].qz > ion_r[2]+HL) vapors[i].qz -= d_size, flag++;
		if (flag>0) {
			vapors[i].px = distvapor(mt) *1e-5;
			vapors[i].py = distvapor(mt) *1e-5;
			vapors[i].pz = distvapor(mt) *1e-5;
		}
	}

}
