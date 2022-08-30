#pragma once
#include "constants.hpp"
#include "flags.hpp"

class Physical{
public:
	double Mion, Mgas, Mvapor, mgas, mvapor, m , m_gas, m_N2, z, myu, alphagas,D0,D0_He,D0_Ar,D0_N2,zzee, cgas, cvapor;
	double c;
	char *gyration_path;
	char *dump_path;
	char *RDF_path;
	double Tnh_gas, Tnh_ion;
	double Ecoeff[3];

// calculation conditions
	int gastype;	/*1:He, 2:Ar, 3:N2*/
	int vaportype;	/*1:MeOH, 2:H2O, 3:EtOH*/
	long int step_relax;
	long int step_repre;
	long int Noftimestep;
	double p;
	double T;
	int Nof_around_gas;
	int Nof_around_vapor;
	int OBSERVE;

	double CG_r0;

	/*******************FUNCTION************************/
	void physicalProp_default(void);
	void PhysicalProp_set(char* condfile, char* file, FLAG *flags);


	PhysicalProp(void);
	~PhysicalProp(void);

private:
	void read(char* infile);
	void set_condition(char* condfile, FLAG *flags, char* file);
	std::vector<double> atype;
};
