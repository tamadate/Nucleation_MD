#pragma once
#include "constants.hpp"

//------------------------------------------------------------------------
class Variables {
public:

	Variables(void);
	~Variables(void){};

	/*variables*/
  int Nth;
	// dammy: CG[0], gas: CG[1], vapor: CG[0]
	std::vector<vector<Molecule>> CG;

	double time;
	double zeta_ion;
	double zeta_gas;

	Potentials Utotal;
	Times times;

	void Uzero(void)	{
		Utotal.Uion=Utotal.Ugas=Utotal.Uvap=Utotal.Ugi=Utotal.Ugg=Utotal.Uvg=Utotal.Uvi=Utotal.Uvv=0;
  }

  void tzero(void)	{times.tion=times.tgas=times.tvap=times.tgi=times.tvv=times.tvg=times.tvi=times.tpair=0;}
	double Usum(void)	{return Utotal.Uion+Utotal.Ugas+Utotal.Uvap+Utotal.Ugi+Utotal.Ugg+Utotal.Uvg+Utotal.Uvi+Utotal.Uvv;}

	/*vectors for potential calculation*/
  std::vector<Pair> ion_pairs;
	std::vector<Pair> pairs_gi;	/*	gas-ion interaction pair list	*/
	std::vector<Pair> pairs_gv;	/*	gas-vapor interaction pair list	*/
	std::vector<Pair> pairs_gg;	/*	gas-gas interaction pair list	*/
	std::vector<Pair> pairs_vv;	/*	gas-gas interaction pair list	*/
	std::vector<Bond> bonds;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;
	std::vector<Atom_type> atypes;
	std::vector<Bond_type> btypes;
	std::vector<Angle_type> ctypes;
	std::vector<Dihedral_type> dtypes;
	std::vector<Atom> atomVapor;
	std::vector<Bond> bonds_v;
	std::vector<Angle> angles_v;
	std::vector<Dihedral> dihedrals_v;
	std::vector<Atom> atomGas(int gastype);
	std::vector<vector<vector<double>>> pair_coeff;

	double bornCoeff[2][2][5];

	void setGasPotentials(void);
	void setBMHPotential(void);
	void setCrossPotentials(int Nion,int Nvapor);

	/*initialization and export to dump file*/
	void read_initial(char* ionFile, char* vaporFile);
	int readIonFile(char* infile);
	int readVaporFile(char* infile);
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
