#include "../md.hpp"

//------------------------------------------------------------//
/*	Make ion-gas interaction pair list  */
//------------------------------------------------------------//

void
MD::make_pairLJCoul(void){
	vars->pairsLJCoul.clear();
	Molecule *mols=vars->Molecules.data();
	for (auto &i : vars->MolID[2]){
		if(vars->Region[i]==CG) continue;
		int j=0;
		for(auto &ion : vars->Molecules[0].inAtoms){
			double dx=mols[i].qx-ion.qx;
			double dy=mols[i].qy-ion.qy;
			double dz=mols[i].qz-ion.qz;
			double r2 = (dx * dx + dy * dy + dz * dz);
			if(r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairsLJ.push_back(p);
			}
			j++;
		}
		i++;
	}
}
