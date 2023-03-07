//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


void
MD::positionLog(void){
	updateVaporinCenters();
	for(auto &i : vars->vapor_in){
		double minDist=1e10;
		int closeAtomID=0;
		for(auto &j : stickPositionList){
			double dx=vars->vapors[i].qx-vars->ions[j].qx;
			double dy=vars->vapors[i].qy-vars->ions[j].qy;
			double dz=vars->vapors[i].qz-vars->ions[j].qz;
			double dr2=dx*dx+dy*dy+dz*dz;
			if(dr2<minDist){
				closeAtomID=j;
				minDist=dr2;
			}
		}
		sprintf(filepath, "stickPositionLog_%d.dat", int(calculation_number));
		FILE*f=fopen(filepath, "a");
		fprintf(f, "%ld\t%d\t%d\t%e\n", itime,i,closeAtomID,sqrt(minDist));
		fclose(f);
	}
}
