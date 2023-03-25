#pragma once
#include "constants.hpp"

//------------------------------------------------------------------------
class Variables {
public:

	int calcID;

	Variables(void);
	~Variables(void){};

	/*variables*/
  	int Nth;
	std::vector<Atom> ions;
	std::vector<Molecule> gases;
  	std::vector<Molecule> vapors;

	double IonX[3];
	double IonV[3];
	double gasX[3];
	double gasV[3];
	double preIonX[3];

	double time;
	double zeta_ion;

	Potentials U;
	Times times;

	void Uzero(void)	{U.Uion=U.Ugas=U.Uvap=U.Ugi=U.Ugg=U.Uvg=U.Uvi=U.Uvv=0;}

  	void tzero(void)	{times.tion=times.tgas=times.tvap=times.tgi=times.tvv=times.tvg=times.tvi=times.tpair=0;}
	double Usum(void)	{return U.Uion+U.Ugas+U.Uvap+U.Ugi+U.Ugg+U.Uvg+U.Uvi+U.Uvv;}

	bool eflag;

	std::vector<int> gas_in;	/*	gas list around ion1	*/
	std::vector<int> gas_out;	/*	gas list far from ion1	*/
	std::vector<int> vapor_in;	/*	gas list around ion1	*/
	std::vector<int> vapor_out;	/*	gas list far from ion1	*/

	/*vectors for potential calculation*/
	std::vector<Bond> bonds;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;
	std::vector<Atom_type> atypes;
	std::vector<Bond_type> btypes;
	std::vector<Angle_type> ctypes;
	std::vector<Dihedral_type> dtypes;
	std::vector<Atom> atomVapor;
	std::vector<Atom> atomGas(int gastype);
	std::vector<vector<vector<double>>> pair_coeff;

	std::vector<Bond> bonds_v;
	std::vector<Angle> angles_v;
	std::vector<Dihedral> dihedrals_v;
	std::vector<Bond_type> btypes_v;
	std::vector<Angle_type> ctypes_v;
	std::vector<Dihedral_type> dtypes_v;

	double bornCoeff[2][2][5];

	void setBMHPotential(void);
	void setCrossPotentials(void);

	/*initialization and export to dump file*/
	void readIonFile(char* infile);
	void readVaporFile(char* infile);
	void ionInitialVelocity(double T);
	void ionRotation(void);
	double totalPotential;
	double totalVirial;

	void ROTATION(double &X, double &Y, double &Z, double A, double B, double C, double x, double y, double z){
		X = cos(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+sin(C)*(y*cos(A)+z*sin(A));
		Y = -sin(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+cos(C)*(y*cos(A)+z*sin(A));
		Z = x*sin(B)-y*sin(A)*cos(B)+z*cos(A)*cos(B);
	}



private:
};
//------------------------------------------------------------------------
