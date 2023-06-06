#pragma once
#include "constants.hpp"

class Physical{
public:
	double Mion, Mgas, Mvapor;	// molecular weights of ion, gas, and vapor
	double mgas, mvapor, m;		// mass of an ion, a gas, and a vapor molecule
	double c;					// mean thermal speed of ion
	double z; 					// number of charges of ion
	double alphagas;			// polarizability of gas molecule
	double zzee;				// q*q
	double T;					// system (gas) temperature
	double p;					// system (gas) pressure

	void setPhysicalProp(int gastype);	// set physical properties
	void readIonProp(char* infile);		// read/set ion properties
	void readVaporProp(char* infile);	// read/set vapor properties

	Physical(void){
		// setting default values
		Mion=Mgas=Mvapor=z=0;
		T=300;
		p=1e5;
	};
	~Physical(void){};

private:

};
