#include "../md.hpp"

/*########################################################################################

-----Initialization-----
intialization:
This intialize the calculation. Reading initial positions and setting initial velocities.
reintialization:
This makes connection between thermal relaxation and main diffusion coeficient calculation. It reset position, time, pair list and margine length.

#######################################################################################*/


void
MD::takeOver(void) {

	ifstream stream(takeOverFile);
	string str;
	int loop=0;
	int Nion=vars->AA[0][0].inAtoms.size();
	int Nvap=vars->atomVapor.size();

	while(getline(stream,str)) {
		if(str.length()==0) continue;
		std::vector<string> readings;
		istringstream stream(str);
		string reading;
		while(getline(stream,reading,'\t')){
			readings.push_back(reading);
		}
		Molecule mol;
		mol.inAtoms=vars->atomVapor;
		for (auto &b : mol.inAtoms){
			for (int thread=0;thread<Nth;thread++){
				b.fxMP.push_back(0);
				b.fyMP.push_back(0);
				b.fzMP.push_back(0);
			}
		}
		mol.bonds=vars->bonds_v;
		mol.angles=vars->angles_v;
		mol.dihedrals=vars->dihedrals_v;

		if(loop<Nion){
			vars->AA[0][0].inAtoms[loop].qx=reading[1];
			vars->AA[0][0].inAtoms[loop].qy=reading[2];
			vars->AA[0][0].inAtoms[loop].qz=reading[3];
			vars->AA[0][0].inAtoms[loop].px=reading[4];
			vars->AA[0][0].inAtoms[loop].py=reading[5];
			vars->AA[0][0].inAtoms[loop].pz=reading[6];
		}
		else{
			int atID=(loop-Nion)%Nvap;
			mol.inAtoms[atID].qx=reading[1];
			mol.inAtoms[atID].qy=reading[2];
			mol.inAtoms[atID].qz=reading[3];
			mol.inAtoms[atID].px=reading[4];
			mol.inAtoms[atID].py=reading[5];
			mol.inAtoms[atID].pz=reading[6];
			if(atID==Nvap-1){
				vars->AA[2].push_back(mol);
			}
		}
		loop++;
	}
	updateInCenters();
}
