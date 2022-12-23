#include "../variables.hpp"

int
Variables::readIonFile(char* infile){
  ifstream stream(infile);
	string str;
	int iflag=0;
	int num_atoms=0;
  int Natypes=atypes.size();
  int Nbtypes=btypes.size();
  int Nctypes=ctypes.size();
  int Ndtypes=dtypes.size();


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

    // Atomic
		if (iflag==1) {
			Atom_type at;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
        if (loop==1) {
          if(tmp=="#c"||tmp=="#cs"||tmp=="#c1"||tmp=="#c2"||tmp=="#c3"||tmp=="#ca"
          ||tmp=="#cp"||tmp=="#cq"||tmp=="#cc"||tmp=="#cd"||tmp=="#ce"||tmp=="#cf"
          ||tmp=="#cg" ||tmp=="#ch"||tmp=="#cx"||tmp=="#cy"||tmp=="#cu"||tmp=="#cv"||tmp=="#cz") at.name="C";
          else if(tmp=="#h1"||tmp=="#h2"||tmp=="#h3"||tmp=="#h4"||tmp=="#h5"||tmp=="#ha"
          ||tmp=="#hc"||tmp=="#hn"||tmp=="#ho"||tmp=="#hp"||tmp=="#hs"||tmp=="#hw"||tmp=="#hx") at.name="H";
          else if(tmp=="f") at.name="F";
          else if(tmp=="cl") at.name="Cl";
          else if(tmp=="br") at.name="Br";
          else if(tmp=="i") at.name="I";
          else if(tmp=="#n"||tmp=="#n1"||tmp=="#n2"||tmp=="#n3"||tmp=="#n4"||tmp=="#na"
          ||tmp=="#nb"||tmp=="#nc"||tmp=="#nd"||tmp=="#ne"||tmp=="#nf"||tmp=="#nh"
          ||tmp=="#no"||tmp=="#ns"||tmp=="#nt"||tmp=="#nx"||tmp=="#ny"||tmp=="#nz"
          ||tmp=="#n+"||tmp=="#nu"||tmp=="#nv"||tmp=="#n7"||tmp=="#n8"||tmp=="#n9") at.name="N";
          else if(tmp=="#o"||tmp=="#oh"||tmp=="#os"||tmp=="#ow") at.name="O";
          else if(tmp=="#p2"||tmp=="#p3"||tmp=="#p4"||tmp=="#p5"||tmp=="#pb"||tmp=="#pc"
          ||tmp=="#pd"||tmp=="#pe"||tmp=="#pf"||tmp=="#px"||tmp=="#py") at.name="P";
          else if(tmp=="#s"||tmp=="#s2"||tmp=="#s4"||tmp=="#s6"||tmp=="#sh"||tmp=="#ss"
          ||tmp=="#sx"||tmp=="#sy") at.name="S";
          else at.name=tmp;
        }
				if (loop==2) at.mass=stod(tmp);
				if (loop==3) at.coeff1=stod(tmp);
				if (loop==4) at.coeff2=stod(tmp);
				loop++;
			}
			atypes.push_back(at);
			num_atoms++;
		}

    // Bond potential parameters
		if (iflag==2) {
			Bond_type bt;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) bt.coeff[0]=stod(tmp);
				if (loop==3) bt.coeff[1]=stod(tmp);
				loop++;
			}
			btypes.push_back(bt);
		}

    // Angle potential parameters
		if (iflag==3) {
			Angle_type at;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) at.coeff[0]=stod(tmp);
				if (loop==3) at.coeff[1]=stod(tmp)/180.0*M_PI;
				loop++;
			}
			ctypes.push_back(at);
		}

    // Dihedral potential parameters
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
			dtypes.push_back(dit);
		}

    // Atom list
		if (iflag==5) {
			Atom a;
			int loop=0;
			a.fx=a.fy=a.fz=a.px=a.py=a.pz=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) a.id=stoi(tmp);
				if (loop==1) {
          a.type=stoi(tmp)-1+Natypes;
          a.mass=atypes[a.type].mass;
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
			Molecules[0].inAtoms.push_back(a);
		}

    // Bond list
		if (iflag==6) {
		 	Bond b;
			int loop=0;
			int atom1, atom2, type;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) b.atom1 = stoi(tmp)-1;
				if (loop==1) b.atom2 = stoi(tmp)-1;
				if (loop==2) b.type=stoi(tmp)-1+Nbtypes;
				loop++;
			}
      Molecules[0].bonds.push_back(b);
		}

    // Angle list
		if (iflag==7) {
			Angle c;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) c.atom1=stoi(tmp)-1;
				if (loop==1) c.atom2=stoi(tmp)-1;
				if (loop==2) c.atom3=stoi(tmp)-1;
				if (loop==3) c.type=stoi(tmp)-1+Nctypes;
				loop++;
			}
			Molecules[0].angles.push_back(c);
		}

    // Dihedral list
		if (iflag==8) {
			Dihedral d;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) d.atom1=stoi(tmp)-1;
				if (loop==1) d.atom2=stoi(tmp)-1;
				if (loop==2) d.atom3=stoi(tmp)-1;
				if (loop==3) d.atom4=stoi(tmp)-1;
				if (loop==4) d.type=stoi(tmp)-1+Ndtypes;
				loop++;
			}
			Molecules[0].dihedrals.push_back(d);
		}
	}

  return num_atoms;
}
