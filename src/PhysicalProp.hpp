#pragma once
#include "constants.hpp"

class Physical{
public:
	double Mion, Mgas, Mvapor;
	double mgas, mvapor, m;
	double c;
	double z, alphagas, zzee;
	double T;
	double p;

	/*******************FUNCTION************************/
	void setPhysicalProp(int gastype);
	void readIonProp(char* infile);
	void readVaporProp(char* infile);

	Physical(void){
		Mion=Mgas=Mvapor=z=0;
		T=300;
		p=1e5;
	};
	~Physical(void){};

private:

};
