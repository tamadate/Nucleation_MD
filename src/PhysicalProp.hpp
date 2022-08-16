#pragma once
#include "CalculationCond.hpp"
#include "flags.hpp"

class Physical{
public:
	const double MN2=28.0, MHe=4.0026, MAr=39.948;
	const double MMeOH=32.0, MH2O=18.0, MEtOH=46.0;
	const double mN2=MN2/Nw/1000.0;
	const double myuHe=1.47e-5;
	const double alphaHe=0.208, alphaN2=1.7*0.5, alphaAr=1.664;
    const double SAr=114, TrefAr=273.0, myuAr=2.125e-5; // REF
    const double SN2=107, TrefN2=273.0, myuN2=1.663e-5; // REF
	double Mion, Mgas, Mvapor, mgas, mvapor, m , m_gas, m_N2, z, myu, alphagas,D0,D0_He,D0_Ar,D0_N2,zzee, cgas, cvapor;
	int num_gas;
	double c;
	char *gyration_path;
	char *dump_path;
	char *RDF_path;
	double Tnh_gas, Tnh_ion;
	int ft;
	double Ecoeff[3];

	long int step_relax;
	long int step_repre;

	/*******************FUNCTION************************/
	void PhysicalProp_set(char* condfile, char* file, FLAG *flags);
	double CG_r0;

private:
	void read(char* infile);
	void set_condition(char* condfile, FLAG *flags, char* file);
	std::vector<double> atype;
};
