#pragma once
#include "constants.hpp"

class Physical{
public:
	double Mion, Mgas, Mvapor;
	double mgas, mvapor, m , m_gas, m_N2;
	double z, myu, alphagas, zzee, cgas, cvapor;
	double c;
	char *gyration_path;
	char *dump_path;
	char *RDF_path;
	double Tnh_gas, Tnh_ion;


	/*******************FUNCTION************************/
	void setPhysicalProp(int gastype, int vaportype, double T, double p);
	void readAtomsFile(char* infile);

	Physical(void){};
	~Physical(void){};

private:
	std::vector<double> atype;
};
