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
	Molecule *mols = vars -> Molecules.data();

	for (auto i : vars->MolID[1]){
		flag=flagx=flagy=flagz=0;
		double dx=mols[i].qx-mols[0].qx;
		double dy=mols[i].qy-mols[0].qy;
		double dz=mols[i].qz-mols[0].qz;
		double dr2=dx*dx+dy*dy+dz*dz;
		mols[i].w=weightFunc(dr2);

		if (dx < -HL) mols[i].qx += d_size, flagx--, flag++;
		if (dy < -HL) mols[i].qy += d_size, flagy--, flag++;
		if (dz < -HL) mols[i].qz += d_size, flagz--, flag++;
		if (dx > HL) mols[i].qx -= d_size, flagx++, flag++;
		if (dy > HL) mols[i].qy -= d_size, flagy++, flag++;
		if (dz > HL) mols[i].qz -= d_size, flagz++, flag++;
		if (flag>0) {
			if(mbdist->number>mbdist->vflux.size()*0.9) {mbdist->makeWeightedMB(pp->cgas,pp->mgas,T);}
			vx = mols[i].px;
	    vy = mols[i].py;
	    vz = mols[i].pz;
	    v2 = vx*vx+vy*vy+vz*vz;
	    v = sqrt(v2);
	    vMB = mbdist->vflux[mbdist->number]*1e-5;
	    mod_factor= vMB/v;
	    mols[i].px = vx * mod_factor;
	    mols[i].py = vy * mod_factor;
	    mols[i].pz = vz * mod_factor;
	    mbdist->number++;
		}
	}
}



void
MD::boundary_scaling_vapor_move(void){
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;
	Molecule *mols = vars -> Molecules.data();

	for (auto i : vars->MolID[2]){
		flag=flagx=flagy=flagz=0;
		double dx=mols[i].qx-mols[0].qx;
		double dy=mols[i].qy-mols[0].qy;
		double dz=mols[i].qz-mols[0].qz;
		double dr2=dx*dx+dy*dy+dz*dz;
		mols[i].w=weightFunc(dr2);

		if (dx < -HL) mols[i].qx += d_size, flagx--, flag++;
		if (dy < -HL) mols[i].qy += d_size, flagy--, flag++;
		if (dz < -HL) mols[i].qz += d_size, flagz--, flag++;
		if (dx > HL) mols[i].qx -= d_size, flagx++, flag++;
		if (dy > HL) mols[i].qy -= d_size, flagy++, flag++;
		if (dz > HL) mols[i].qz -= d_size, flagz++, flag++;
		if (flag>0) {
			if(mbdistV->number>mbdistV->vflux.size()*0.9) {mbdistV->makeWeightedMB(pp->cvapor,pp->mvapor,T);}
	    vx=mols[i].px;
	    vy=mols[i].py;
	    vz=mols[i].pz;
			v2 = vx*vx+vy*vy+vz*vz;
			v = sqrt(v2);
			vMB = mbdistV->vflux[mbdistV->number]*1e-5;
			mod_factor= vMB/v;
	    mols[i].px = vx * mod_factor;
	    mols[i].py = vy * mod_factor;
	    mols[i].pz = vz * mod_factor;
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
	Molecule *mols = vars -> Molecules.data();

	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	for (auto i : vars->MolID[1]){
		flag=flagx=flagy=flagz=0;
		double dx=mols[i].qx-mols[0].qx;
		double dy=mols[i].qy-mols[0].qy;
		double dz=mols[i].qz-mols[0].qz;
		if (dx < -HL) mols[i].qx += d_size, flagx--, flag++;
		if (dy < -HL) mols[i].qy += d_size, flagy--, flag++;
		if (dz < -HL) mols[i].qz += d_size, flagz--, flag++;
		if (dx > HL) mols[i].qx -= d_size, flagx++, flag++;
		if (dy > HL) mols[i].qy -= d_size, flagy++, flag++;
		if (dz > HL) mols[i].qz -= d_size, flagz++, flag++;
		if (flag>0) {
			mols[i].px = distgas(mt) *1e-5;
			mols[i].py = distgas(mt) *1e-5;
			mols[i].pz = distgas(mt) *1e-5;
		}
	}

	normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));
	for (auto i : vars->MolID[2]){
		flag=flagx=flagy=flagz=0;
		double dx=mols[i].qx-mols[0].qx;
		double dy=mols[i].qy-mols[0].qy;
		double dz=mols[i].qz-mols[0].qz;
		if (mols[i].qx < mols[0].qx-HL) mols[i].qx += d_size, flag++;
		if (mols[i].qy < mols[0].qy-HL) mols[i].qy += d_size, flag++;
		if (mols[i].qz < mols[0].qz-HL) mols[i].qz += d_size, flag++;
		if (mols[i].qx > mols[0].qx+HL) mols[i].qx -= d_size, flag++;
		if (mols[i].qy > mols[0].qy+HL) mols[i].qy -= d_size, flag++;
		if (mols[i].qz > mols[0].qz+HL) mols[i].qz -= d_size, flag++;
		if (flag>0) {
			mols[i].px = distvapor(mt) *1e-5;
			mols[i].py = distvapor(mt) *1e-5;
			mols[i].pz = distvapor(mt) *1e-5;
		}
	}

}
