#include "../md.hpp"
#include <random>
#include <algorithm>
#include <dirent.h>
#include <vector>



void
MD::boundary_scaling_gas_move(void){
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;

	for (auto &a : vars->CG[1]){
		flag=flagx=flagy=flagz=0;
		if (a.qx < vars->AA[0][0].qx-HL) a.qx += d_size, flagx--, flag++;
		if (a.qy < vars->AA[0][0].qy-HL) a.qy += d_size, flagy--, flag++;
		if (a.qz < vars->AA[0][0].qz-HL) a.qz += d_size, flagz--, flag++;
		if (a.qx > vars->AA[0][0].qx+HL) a.qx -= d_size, flagx++, flag++;
		if (a.qy > vars->AA[0][0].qy+HL) a.qy -= d_size, flagy++, flag++;
		if (a.qz > vars->AA[0][0].qz+HL) a.qz -= d_size, flagz++, flag++;
		if (flag>0) {
			if(mbdist->number>mbdist->vflux.size()*0.9) {mbdist->makeWeightedMB(pp->cgas,pp->mgas,T);}
			vx = a.px;
	    vy = a.py;
	    vz = a.pz;
	    v2 = vx*vx+vy*vy+vz*vz;
	    v = sqrt(v2);
	    vMB = mbdist->vflux[mbdist->number]*1e-5;
	    mod_factor= vMB/v;
	    a.px = vx * mod_factor;
	    a.py = vy * mod_factor;
	    a.pz = vz * mod_factor;
	    mbdist->number++;
		}
	}
}



void
MD::boundary_scaling_vapor_move(void){
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;

	for (auto &a : vars->CG[2]){
		flag=flagx=flagy=flagz=0;
		if (a.qx < vars->AA[0][0].qx-HL) a.qx += d_size, flagx--, flag++;
		if (a.qy < vars->AA[0][0].qy-HL) a.qy += d_size, flagy--, flag++;
		if (a.qz < vars->AA[0][0].qz-HL) a.qz += d_size, flagz--, flag++;
		if (a.qx > vars->AA[0][0].qx+HL) a.qx -= d_size, flagx++, flag++;
		if (a.qy > vars->AA[0][0].qy+HL) a.qy -= d_size, flagy++, flag++;
		if (a.qz > vars->AA[0][0].qz+HL) a.qz -= d_size, flagz++, flag++;
		if (flag>0) {
			if(mbdistV->number>mbdistV->vflux.size()*0.9) {mbdistV->makeWeightedMB(pp->cvapor,pp->mvapor,T);}
	    vx=a.px;
	    vy=a.py;
	    vz=a.pz;
			v2 = vx*vx+vy*vy+vz*vz;
			v = sqrt(v2);
			vMB = mbdistV->vflux[mbdistV->number]*1e-5;
			mod_factor= vMB/v;
	    a.px = vx * mod_factor;
	    a.py = vy * mod_factor;
	    a.pz = vz * mod_factor;
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

	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	for (auto &a : vars->CG[1]){
		flag=flagx=flagy=flagz=0;
		if (a.qx < vars->AA[0][0].qx-HL) a.qx += d_size, flagx--, flag++;
		if (a.qy < vars->AA[0][0].qy-HL) a.qy += d_size, flagy--, flag++;
		if (a.qz < vars->AA[0][0].qz-HL) a.qz += d_size, flagz--, flag++;
		if (a.qx > vars->AA[0][0].qx+HL) a.qx -= d_size, flagx++, flag++;
		if (a.qy > vars->AA[0][0].qy+HL) a.qy -= d_size, flagy++, flag++;
		if (a.qz > vars->AA[0][0].qz+HL) a.qz -= d_size, flagz++, flag++;
		if (flag>0) {
			a.px = distgas(mt) *1e-5;
			a.py = distgas(mt) *1e-5;
			a.pz = distgas(mt) *1e-5;
		}
	}

	normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));
	for (auto &a : vars->CG[2]){
		flag=flagx=flagy=flagz=0;
		if (a.qx < vars->AA[0][0].qx-HL) a.qx += d_size, flag++;
		if (a.qy < vars->AA[0][0].qy-HL) a.qy += d_size, flag++;
		if (a.qz < vars->AA[0][0].qz-HL) a.qz += d_size, flag++;
		if (a.qx > vars->AA[0][0].qx+HL) a.qx -= d_size, flag++;
		if (a.qy > vars->AA[0][0].qy+HL) a.qy -= d_size, flag++;
		if (a.qz > vars->AA[0][0].qz+HL) a.qz -= d_size, flag++;
		if (flag>0) {
			a.px = distvapor(mt) *1e-5;
			a.py = distvapor(mt) *1e-5;
			a.pz = distvapor(mt) *1e-5;
		}
	}

}
