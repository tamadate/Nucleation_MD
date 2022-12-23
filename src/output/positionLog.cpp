//------------------------------------------------------------------------
#include "../md.hpp"
//------------------------------------------------------------------------


void
MD::positionLog(void){
	updateInCenters();
 	Molecule *mols=vars->Molecules.data();
	Atom *ion=vars->Molecules[0].inAtoms.data();
	int i=0;
	for(auto i : vars->MolID[2]){
		double minDist=1e10;
		int closeAtomID=0;
		for(auto &j : stickPositionList){
			double dx=mols[i].qx-ion[j].qx;
			double dy=mols[i].qy-ion[j].qy;
			double dz=mols[i].qz-ion[j].qz;
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
		i++;
	}
}
