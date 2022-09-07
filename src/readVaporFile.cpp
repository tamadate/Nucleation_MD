#include "variables.hpp"

int
Variables::readVaporFile(char* infile){

  ifstream stream(infile);
	string str;
	int iflag=0;
	int num_atoms=0;

	while(getline(stream,str)) {
		if(str.length()==0) continue;
		if (str=="atom type name mass coeff1 coeff2") {iflag=1; continue;}
		if (str=="bond type name coeff1 coeff2") {iflag=2; continue;}
		if (str=="angle type name coeff1 coeff2") {iflag=3; continue;}
		if (str=="dihedral type name coeff1 coeff2 coeff3 coeff4") {iflag=4; continue;}
		if (str=="atoms") {iflag=5; continue;}
		if (str=="bonds") {iflag=6; continue;}
		if (str=="angles") {iflag=7; continue;}
		if (str=="dihedrals") {iflag=8; continue;}
	    string tmp;
    	istringstream stream(str);
		if (iflag==1) {
			Atom_type at;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) at.mass=stod(tmp);
				if (loop==3) at.coeff1=stod(tmp);
				if (loop==4) at.coeff2=stod(tmp);
				loop++;
			}
			atypes_v.push_back(at);
			num_atoms++;
		}
		if (iflag==2) {
			Bond_type bt;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) bt.coeff[0]=stod(tmp);
				if (loop==3) bt.coeff[1]=stod(tmp);
				loop++;
			}
			btypes_v.push_back(bt);
		}
		if (iflag==3) {
			Angle_type at;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) at.coeff[0]=stod(tmp);
				if (loop==3) at.coeff[1]=stod(tmp)/180.0*M_PI;
				loop++;
			}
			ctypes_v.push_back(at);
		}
		if (iflag==4) {
			Dihedral_type dit;
			int loop=0;
			double coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7, coeff8, coeff9, coeff10,  coeff11, coeff12, coeff13, coeff14, coeff15;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) dit.multi=stoi(tmp);
				if (loop==3) dit.coeff[0]=stod(tmp);
				if (loop==4) dit.coeff[1]=stod(tmp);
				if (loop==5) {
					dit.coeff[2]=stod(tmp);
					dit.coeff[3]=cos(stod(tmp)/180.0*M_PI);
					dit.coeff[4]=sin(stod(tmp)/180.0*M_PI);
				}
				if (loop==6) dit.coeff[5]=stod(tmp);
				if (loop==7) dit.coeff[6]=stod(tmp);
				if (loop==8) {
					dit.coeff[7]=stod(tmp);
					dit.coeff[8]=cos(stod(tmp)/180.0*M_PI);
					dit.coeff[9]=sin(stod(tmp)/180.0*M_PI);
				}
				if (loop==9) dit.coeff[10]=stod(tmp);
				if (loop==10) dit.coeff[11]=stod(tmp);
				if (loop==11) {
					dit.coeff[12]=stod(tmp);
					dit.coeff[13]=cos(stod(tmp)/180.0*M_PI);
					dit.coeff[14]=sin(stod(tmp)/180.0*M_PI);
				}
				loop++;
			}
			dtypes_v.push_back(dit);
		}
		if (iflag==5) {
			Atom a;
			int loop=0;
			a.fx=a.fy=a.fz=a.px=a.py=a.pz=0;
			for (int thread=0;thread<Nth;thread++){
				a.fxMP.push_back(0);
				a.fyMP.push_back(0);
				a.fzMP.push_back(0);
			}
			while(getline(stream,tmp,'\t')) {
				if (loop==0) a.id=stoi(tmp);
				if (loop==1) {
          a.type=stoi(tmp)-1;
          a.mass=atypes_v[a.type].mass;
        }
				if (loop==2) a.charge=stod(tmp);
				if (loop==3) a.qx=stod(tmp);
				if (loop==4) a.qy=stod(tmp);
				if (loop==5) a.qz=stod(tmp);
				if (loop==6) a.px=stod(tmp);
				if (loop==7) a.py=stod(tmp);
				if (loop==8) a.pz=stod(tmp);
				loop++;
			}
			atomVapor.push_back(a);
		}
		if (iflag==6) {
		 	Bond b;
			int loop=0;
			int atom1, atom2, type;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) b.atom1 = stoi(tmp)-1;
				if (loop==1) b.atom2 = stoi(tmp)-1;
				if (loop==2) b.type=stoi(tmp)-1;
				loop++;
			}
			bonds_v.push_back(b);
		}
		if (iflag==7) {
			Angle c;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) c.atom1=stoi(tmp)-1;
				if (loop==1) c.atom2=stoi(tmp)-1;
				if (loop==2) c.atom3=stoi(tmp)-1;
				if (loop==3) c.type=stoi(tmp)-1;
				loop++;
			}
			angles_v.push_back(c);
		}
		if (iflag==8) {
			Dihedral d;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) d.atom1=stoi(tmp)-1;
				if (loop==1) d.atom2=stoi(tmp)-1;
				if (loop==2) d.atom3=stoi(tmp)-1;
				if (loop==3) d.atom4=stoi(tmp)-1;
				if (loop==4) d.type=stoi(tmp)-1;
				loop++;
			}
			dihedrals_v.push_back(d);
		}
	}
	pair_coeff.resize(num_atoms);
	for (int i=0;i<num_atoms;i++){
		pair_coeff_v[i].resize(num_atoms);
		for (int j=0;j<num_atoms;j++){
			pair_coeff_v[i][j].resize(2);
		}
	}
	for (int i=0;i<num_atoms;i++){
		for (int j=0;j<num_atoms;j++){
			double epu=sqrt(atypes_v[i].coeff1*atypes_v[j].coeff1);
			double sigma=(atypes_v[i].coeff2+atypes_v[j].coeff2)*0.5;
			//sigma=sqrt(atypes_v[i].coeff2*atypes_v[j].coeff2);
			pair_coeff_v[i][j][0]=48 * epu*pow(sigma,12.0);
			pair_coeff_v[j][j][1]=24 * epu*pow(sigma,6.0);
		}
	}
  return num_atoms;
}
