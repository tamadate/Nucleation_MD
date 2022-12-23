#include "../md.hpp"


//------------------------------------------------------------//
/*	Make update regions  */
//------------------------------------------------------------//
void
MD::updateInOut(void){
	Molecule *mols=vars->Molecules.data();
	for (auto i : vars->MolID[1]){
		double dx = mols[i].qx - mols[0].qx;
		double dy = mols[i].qy - mols[0].qy;
		double dz = mols[i].qz - mols[0].qz;
		double r2 = (dx * dx + dy * dy + dz * dz);
		double dr2=dx*dx+dy*dy+dz*dz;
		int oriFlag=vars->Region[i];
		if(dr2<RO2) vars->Region[i]=AA;
		else {
			if(dr2<RI2) {
				vars->Region[i]=AACG;
				if(oriFlag==CG) makeDiatomicProp_in(mols[i]);
			}
			else vars->Region[i]=CG;
		}
	}
	for (auto i : vars->MolID[2]){
		double dx = mols[i].qx - mols[0].qx;
		double dy = mols[i].qy - mols[0].qy;
		double dz = mols[i].qz - mols[0].qz;
		double r2 = (dx * dx + dy * dy + dz * dz);
		double dr2=dx*dx+dy*dy+dz*dz;
		int oriFlag=vars->Region[i];
		if(dr2<RO2) vars->Region[i]=AA;
		else {
			if(dr2<RI2) {
				vars->Region[i]=AACG;
				if(oriFlag==CG) makePolyatomicProp_in(mols[i]);
			}
			else {
				vars->Region[i]=CG;
				if(oriFlag!=CG){
					sprintf(filepath, "vapor_out_%d.dat", int(calculation_number));
					FILE*f=fopen(filepath, "a");
					fprintf(f, "%d %e %e %e %e %e %e %e\n", mols[0].id,itime*dt,dx,dy,dz,\
					mols[i].px-mols[0].px,mols[i].py-mols[0].py,mols[i].pz-mols[0].pz);
					fclose(f);
				}
			}
		}
	}
}
