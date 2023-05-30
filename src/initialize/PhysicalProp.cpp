#include "../PhysicalProp.hpp"

void
Physical::setPhysicalProp(int gastype){
	if(gastype==1) {
		Mgas=MHe;
		alphagas=alphaHe;
	}
	if(gastype==2) {
		Mgas=MN2;
		alphagas=alphaN2/2.0;
	}
	if(gastype==3) {
		Mgas=MN2;
		alphagas=alphaN2;
	}
	if(gastype==4) {
		Mgas=MAr;
		alphagas=alphaAr;
	}
	mvapor=Mvapor/Nw/1000.0; // mass of a vapor molecule
	mgas=Mgas/Nw/1000.0;	// mass of a gas molecule
	m=Mion/Nw/1000.0;		// mass of an ion
	c=sqrt(8*kb*T/M_PI/m);	// mean thermal speed of ion

	printf("Ion mass\t\t%f g/mol\nIon charges\t\t%f\n", Mion,z);
	printf("Vapor mass\t\t%f g/mol\n", Mvapor);
	cout<<"**************************************************"<<endl;
	cout<<"**************************************************"<<endl;
	cout<<"**************************************************"<<endl;
}


void
Physical::readIonProp(char* infile){
	ifstream stream(infile);
	string str;
	int iflag=0;
	Mion=0;
	std::vector<double> mass_array;
	while(getline(stream,str)) {
		if(str.length()==0) continue;
		if (str=="atom type") {iflag=1; continue;}
		if (str=="atoms") {iflag=2; continue;}
		if (str=="bonds") {break;}
	    string tmp;
    	istringstream stream(str);
		if (iflag==1) {
			int loop=0;
			double mass;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) mass=stod(tmp);
				loop++;
			}
			mass_array.push_back(mass);
		}
		if (iflag==2) {
			int loop=0;
			int type, id;
			double charge;
			while(getline(stream,tmp,'\t')) {
				if (loop==1) type=stoi(tmp);
				if (loop==2) charge=stod(tmp);
				loop++;
			}
			Mion+=mass_array[type-1];
			z+=charge;
		}
	}
}

void
Physical::readVaporProp(char* infile){
	ifstream stream(infile);
	string str;
	int iflag=0;
	Mvapor=0;
	std::vector<double> atype;
	std::vector<double> atype_v;
	while(getline(stream,str)) {
		if(str.length()==0) continue;
		if (str=="atom type") {iflag=1; continue;}
		if (str=="atoms") {iflag=2; continue;}
		if (str=="bonds") {break;}
	    string tmp;
    	istringstream stream(str);
		if (iflag==1) {
			int loop=0;
			double mass;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) mass=stod(tmp);
				loop++;
			}
			atype_v.push_back(mass);
		}
		if (iflag==2) {
			int loop=0;
			int type, id;
			while(getline(stream,tmp,'\t')) {
				if (loop==1) type=stoi(tmp)-1;
				loop++;
			}
			double mass=atype_v[type];
			Mvapor+=mass;
		}
	}
}
