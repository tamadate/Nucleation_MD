#pragma once
#include "constants.hpp"

//------------------------------------------------------------------------
class Variables {
public:

	int calcID;
	Times times;	// Timer

	Variables(void);
	~Variables(void){};

	/*variables*/
	std::vector<Atom> ions;	// ion atoms array
	std::vector<Molecule> gases;	// gas molecules array
  	std::vector<Molecule> vapors;	// vapor molecules array

	std::vector<int> gas_in;	//	gas list around ion	
	std::vector<int> gas_out;	//	gas list far from ion
	std::vector<int> vapor_in;	//	gas list around ion
	std::vector<int> vapor_out;	//	gas list far from ion

	double IonX[3];	// ion center of mass cordinate
	double IonV[3];	// ion center of mass velocity
	double gasX[3]; // gas center of mass cordinate
	double gasV[3]; // gas center of mass velocity
	double preIonX[3]; // ion center of mass cordinate at previous time step

	double time;

	Potentials U;
	bool eflag;	// flag for potential calculation
	void Uzero(void)	{U.Uion=U.Ugas=U.Uvap=U.Ugi=U.Ugg=U.Uvg=U.Uvi=U.Uvv=0;} // set all potential to zero
	double Usum(void)	{return U.Uion+U.Ugas+U.Uvap+U.Ugi+U.Ugg+U.Uvg+U.Uvi+U.Uvv;} // compute total potential
	double totalPotential;
	double totalVirial;

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

	void setCrossPotentials(void);

	/*initialization and export to dump file*/
	void readIonFile(char* infile);
	void readVaporFile(char* infile);
	void ionInitialVelocity(double T);
	void ionRotation(void);

	void ROTATION(double &X, double &Y, double &Z, double A, double B, double C, double x, double y, double z){
		X = cos(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+sin(C)*(y*cos(A)+z*sin(A));
		Y = -sin(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+cos(C)*(y*cos(A)+z*sin(A));
		Z = x*sin(B)-y*sin(A)*cos(B)+z*cos(A)*cos(B);
	}



private:
};
//------------------------------------------------------------------------
