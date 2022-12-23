#include "../md.hpp"

//------------------------------------------------------------//
/*	Make vapor-vapor interaction pair list */
//------------------------------------------------------------//
void
MD::make_pairsLJCoulHybrid(void){
	vars->pairsLJCoulHybrid.clear();
	Molecule *mols=vars->Molecules.data();

	// make gas-gas pair list
	int vs=vars->MolID[2].size();
	for (int i=0; i<vs-1; i++){
		for (int j=i+1; j<vs; j++){
			double dx = mols[i].qx - mols[j].qx;
			double dy = mols[i].qy - mols[j].qy;
			double dz = mols[i].qz - mols[j].qz;
			adjust_periodic(dx, dy, dz, d_size);
			double r2 = (dx * dx + dy * dy + dz * dz);
			// 00000001|00000001=00000001
			// 00000011|00000001=00000011
			// 00000011|00000011=00000011
			// 00000010|00000011=00000011
			// 00000010|00000010=00000010 ... only this case, no AA interaction

			if((vars->Region[i]|vars->Region[j])==CG){	// if both of two are CG region
				if (r2 < ML2){
					Pair p;
					p.i=i;
					p.j=j;
					vars->pairsLJHybrid.push_back(p);
				}
			}
			else{
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairsLJCoulHybrid.push_back(p);
			}
		}
	}
}
