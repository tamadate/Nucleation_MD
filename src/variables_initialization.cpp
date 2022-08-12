//------------------------------------------------------------------------
#include "variables.hpp"
//------------------------------------------------------------------------
void 
Variables::read_initial(char* infile) {
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
		if (str=="molecules") {iflag=9; continue;}
	    string tmp;
    	istringstream stream(str);		
		if (iflag==1) {
			int loop=0;
			int type;
			double mass, coeff1, coeff2;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) type=stoi(tmp);
				if (loop==2) mass=stod(tmp);
				if (loop==3) coeff1=stod(tmp);
				if (loop==4) coeff2=stod(tmp);
				loop++;
			}
			add_atype(type, mass, coeff1, coeff2);
			num_atoms++;
		}
		if (iflag==2) {
			Bond_type bt;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) bt.type=stoi(tmp);
				if (loop==2) bt.coeff[0]=stod(tmp);
				if (loop==3) bt.coeff[1]=stod(tmp);
				loop++;
			}
			btypes.push_back(bt);
		}
		if (iflag==3) {
			Angle_type at;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) at.type=stoi(tmp);
				if (loop==2) at.coeff[0]=stod(tmp);
				if (loop==3) at.coeff[1]=stod(tmp)/180.0*M_PI;
				loop++;
			}
			ctypes.push_back(at);
		}
		if (iflag==4) {
			Dihedral_type dit;
			int loop=0;
			double coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7, coeff8, coeff9, coeff10,  coeff11, coeff12, coeff13, coeff14, coeff15;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) dit.type=stoi(tmp);
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
		if (iflag==5) {
			int loop=0;
			int type, id;
			double charge, x, y, z, vx, vy, vz, mass;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) id=stoi(tmp);
				if (loop==1) type=stoi(tmp);
				if (loop==2) charge=stod(tmp);
				if (loop==3) x=stod(tmp);
				if (loop==4) y=stod(tmp);
				if (loop==5) z=stod(tmp);
				if (loop==6) vx=stod(tmp);
				if (loop==7) vy=stod(tmp);
				if (loop==8) vz=stod(tmp);
				loop++;
			}
			for (auto &at : atypes) {
				if(at.type==type) mass=at.mass;
			}
			add_ions(id,type,x,y,z,vx,vy,vz,0,0,0,charge,mass);
		}
		if (iflag==6) {
			int loop=0;
			int atom1, atom2, type;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) atom1=stoi(tmp);
				if (loop==1) atom2=stoi(tmp);
				if (loop==2) type=stoi(tmp);
				loop++;
			}
			add_bonds(atom1-1,atom2-1,type-1);
		}
		if (iflag==7) {
			int loop=0;
			int atom1, atom2, atom3, type;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) atom1=stoi(tmp);
				if (loop==1) atom2=stoi(tmp);
				if (loop==2) atom3=stoi(tmp);
				if (loop==3) type=stoi(tmp);
				loop++;
			}
			add_angles(atom1-1,atom2-1,atom3-1,type-1);
		}
		if (iflag==8) {
			int loop=0;
			int atom1, atom2, atom3, atom4, type;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) atom1=stoi(tmp);
				if (loop==1) atom2=stoi(tmp);
				if (loop==2) atom3=stoi(tmp);
				if (loop==3) atom4=stoi(tmp);
				if (loop==4) type=stoi(tmp);
				loop++;
			}
			add_dihedrals(atom1-1,atom2-1,atom3-1,atom4-1,type-1);
		}
		if (iflag==9) {
			int atom;
			atom=stoi(str);
			molecules.push_back(atom);
		}
	}
	pair_coeff.resize(num_atoms+1);
	for (int i=0;i<num_atoms+1;i++){
		pair_coeff[i].resize(num_atoms+1);
		for (int j=0;j<num_atoms+1;j++){
			pair_coeff[i][j].resize(2);
		}
	}
	for (int i=0;i<num_atoms;i++){
		for (int j=0;j<num_atoms;j++){
			double epu, sigma;
			epu=sqrt(atypes[i].coeff1*atypes[j].coeff1);
			//sigma=(atypes[i].coeff2+atypes[j].coeff2)*0.5;
			sigma=sqrt(atypes[i].coeff2*atypes[j].coeff2);
			pair_coeff[atypes[i].type][atypes[j].type][0]=48 * epu*pow(sigma,12.0);
			pair_coeff[atypes[i].type][atypes[j].type][1]=24 * epu*pow(sigma,6.0);
		}
	}

//difine Born-Mayer-Huggins conefficients for "only" NaCl
//https://doi.org/10.1063/1.1522375
//https://doi.org/10.1016/0022-3697(64)90160-X
//0:A/rho, 1:6C, 2:8D, 3:sigma, 4:1/rho
//NaNa, A=25.4435kJ/mol, C=101.1719kJ/mol, D=48.1771kJ/mol, sigma=2.340A, 1/rho=3.1546A-1
//NaCl, A=20.3548kJ/mol, C=674.4793kJ/mol, D=837.077kJ/mol, sigma=2.755A, 1/rho=3.1546A-1
//ClCl, A=15.2661kJ/mol, C=6985.6786kJ/mol, D=14031.5785kJ/mol, sigma=3.170A, 1/rho=3.1546A-1
	bornCoeff[0][0][0]=25.4435*3.1546/4.184;
	bornCoeff[0][1][0]=bornCoeff[1][0][0]=20.3548*3.1546/4.184;
	bornCoeff[1][1][0]=15.2661*3.1546/4.184;

	bornCoeff[0][0][1]=6*101.1719/4.184;
	bornCoeff[0][1][1]=bornCoeff[1][0][1]=6*674.4793/4.184;
	bornCoeff[1][1][1]=6*6985.6786/4.184;

	bornCoeff[0][0][2]=8*48.1771/4.184;
	bornCoeff[0][1][2]=bornCoeff[1][0][2]=8*837.077/4.184;
	bornCoeff[1][1][2]=8*14031.5785/4.184;

	bornCoeff[0][0][3]=2.340;
	bornCoeff[0][1][3]=bornCoeff[1][0][3]=2.755;
	bornCoeff[1][1][3]=3.170;

	bornCoeff[0][0][4]=3.1546;
	bornCoeff[0][1][4]=bornCoeff[1][0][4]=3.1546;
	bornCoeff[1][1][4]=3.1546;


	random_device seed;
	double A,B,C,x,y,z;
	A=seed(),B=seed(),C=seed()
;
	for(auto &a : ions) {x=a.qx,y=a.qy,z=a.qz; ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z);}
    int is=ions.size();
    for(int i=0;i<is-1;i++) {
        for(int j=i+1;j<is;j++){
            int flag=0;
            for (auto &d : dihedrals) {
                int I=d.atom1;
                int J=d.atom2;
                int K=d.atom3;
                int L=d.atom4;
                if(i==I){if(j==J||j==K||j==L) flag=1;}
                if(i==J){if(j==I||j==K||j==L) flag=1;}
                if(i==K){if(j==J||j==I||j==L) flag=1;}
                if(i==L){if(j==J||j==K||j==I) flag=1;}
            }
            if (flag==0){
                Pair p;
                p.i=i;
                p.j=j;
                ion_pairs.push_back(p);
            }
        }
    }
    
}

//------------------------------------------------------------------------
void
Variables::set_initial_velocity(Physical *pp) {
	random_device seed;
	default_random_engine engine(seed());
    for(auto &a : ions) {
        double matom=a.mass/6.02e23/1000.0;
        normal_distribution<> dist(0.0, sqrt(kb*T/matom));
		a.px=dist(engine)*1e-5;
		a.py=dist(engine)*1e-5;
		a.pz=dist(engine)*1e-5;
	}
}
