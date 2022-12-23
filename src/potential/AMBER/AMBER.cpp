#include "potentialAMBER.hpp"
/*########################################################################################

compute intramolecular interaction - AMBER

#######################################################################################*/

void
PotentialAMBER::compute(Variables *vars, FLAG *flags) {
	vars->times.tion-=omp_get_wtime();
	computeLong(vars,flags);
	computeBond(vars,flags);
	computeAngle(vars,flags);
	computeDihedral(vars,flags);
	vars->times.tion+=omp_get_wtime();
}

void
PotentialAMBER::initialAMBER(Variables *vars, FLAG *flags){
	std::vector<Pair> noLong;
	Pair p;
	for (auto &b : vars-> Molecules[0].bonds) {
		p.i=b.atom1;
		p.j=b.atom2;
		if(p.j<p.i) {
			p.i=b.atom2;
			p.j=b.atom1;
		}
		noLong.push_back(p);
	}
	for (auto &b : vars-> Molecules[0].angles) {
		p.i=b.atom1;
		p.j=b.atom3;
		if(p.j<p.i) {
			p.i=b.atom3;
			p.j=b.atom1;
		}
		noLong.push_back(p);
	}
	for (auto &b : vars-> Molecules[0].dihedrals) {
		p.i=b.atom1;
		p.j=b.atom4;
		if(p.j<p.i) {
			p.i=b.atom4;
			p.j=b.atom1;
		}
		noLong.push_back(p);
	}

	Atom *ions=vars->Molecules[0].inAtoms.data();
	const int is=vars->Molecules[0].inAtoms.size();
	for(int i=0;i<is-1;i++){
		for(int j=i+1;j<is;j++){
			int flag=1;
			for(auto &a : noLong){
				if(a.i==i && a.j==j){
					flag=0;
					break;
				}
			}
			if(flag==1) {
				p.i=i;
				p.j=j;
				longPair.push_back(p);
			}
		}
	}
}
