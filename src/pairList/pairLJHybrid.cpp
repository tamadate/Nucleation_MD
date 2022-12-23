
#include "../md.hpp"

//------------------------------------------------------------//
/*	Make gas-gas & gas-vapor interactions pair list */
//------------------------------------------------------------//
void
MD::make_pairLJHybrid(void){
	vars->pairsLJHybrid.clear();
	Molecule *mols=vars->Molecules.data();

	// make ion-gas pair list
	for (auto i : vars->MolID[1]){
		for (auto j : vars->MolID[2]){
			double dx = mols[i].qx - mols[j].qx;
			double dy = mols[i].qy - mols[j].qy;
			double dz = mols[i].qz - mols[j].qz;
			adjust_periodic(dx, dy, dz, d_size);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairsLJHybrid.push_back(p);
			}
		}
	}

	// make gas-gas pair list
	int gs=vars->MolID[1].size();
	for (int i=0; i<gs-1; i++){
		for (int j=i+1; j<gs; j++){
			double dx = mols[i].qx - mols[j].qx;
			double dy = mols[i].qy - mols[j].qy;
			double dz = mols[i].qz - mols[j].qz;
			adjust_periodic(dx, dy, dz, d_size);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairsLJHybrid.push_back(p);
			}
		}
	}
}
