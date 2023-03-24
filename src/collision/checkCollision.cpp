#include "collision.hpp"

void
Collision::checkCollision(int itime){
    for (auto &i : vars->vapor_in) {
		double dmin=10000000;
		for (auto &b : vars->vapors[i].inAtoms){
			for (auto &a : vars->ions) {
				double dx = b.qx - a.qx;
				double dy = b.qy - a.qy;
				double dz = b.qz - a.qz;
				double r2 = (dx * dx + dy * dy + dz * dz);
				double r = sqrt(r2);
				double sig=vars->pair_coeff_vi_se[b.type][a.type][0];
				double d=r-sig;
				if(dmin>d) dmin=d;
			}
		}
		if (dmin > 10 && collisionFlagVapor[i]!=0) {
			outputOut(i);
			collisionFlagVapor[i]=0;
		}
		if (dmin < 0 && collisionFlagVapor[i]==0){
			outputIn(i);
			collisionFlagVapor[i]=itime;
		}
	}
}