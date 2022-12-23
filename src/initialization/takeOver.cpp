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
	int Nion=vars->Molecules[0].inAtoms.size();
	int Nvap=vars->atomVapor.size();
	int Nsofar=vars->Molecules.size();

	int i=0;
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
		mol.bonds=vars->bonds_v;
		mol.angles=vars->angles_v;
		mol.dihedrals=vars->dihedrals_v;

		if(loop<Nion){
			vars->Molecules[0].inAtoms[loop].qx=reading[1];
			vars->Molecules[0].inAtoms[loop].qy=reading[2];
			vars->Molecules[0].inAtoms[loop].qz=reading[3];
			vars->Molecules[0].inAtoms[loop].px=reading[4];
			vars->Molecules[0].inAtoms[loop].py=reading[5];
			vars->Molecules[0].inAtoms[loop].pz=reading[6];
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
				vars->Molecules.push_back(mol);
				vars->MolID[2].push_back(i+Nsofar);
				i++;
			}
		}
		loop++;
	}
	updateInCenters();
}
