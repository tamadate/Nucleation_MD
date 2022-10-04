#pragma once
#include "constants.hpp"

//------------------------------------------------------------------------
class Variables {
public:

	Variables(void);
	~Variables(void){};

	/*variables*/
  int Nth;
	// ion: effectiveIn[0], gas: effectiveIn[1], vapor: effectiveIn[2]
	std::vector<vector<Molecule>> effectiveIn;
	// dammy: effectiveOut[0], gas: effectiveOut[1], vapor: effectiveOut[0]
	std::vector<vector<Molecule_out>> effectiveOut;

	double time;
	double zeta_ion;
	double zeta_gas;

	Potentials Utotal;
  std::vector<Potentials> U_MP;
	Times times;

	void Uzero(void)	{
		Utotal.Uion=Utotal.Ugas=Utotal.Uvap=Utotal.Ugi=Utotal.Ugg=Utotal.Uvg=Utotal.Uvi=Utotal.Uvv=0;
		for (int nth=0;nth<Nth;nth++){
			U_MP[nth].Uion=U_MP[nth].Ugas=U_MP[nth].Uvap=U_MP[nth].Ugi=U_MP[nth].Ugg=U_MP[nth].Uvg=U_MP[nth].Uvi=U_MP[nth].Uvv=0;
		}
  }
  void Ucombine(void)	{
    for (int nth=0;nth<Nth;nth++){
      Utotal.Uion+=U_MP[nth].Uion;
      Utotal.Uvap+=U_MP[nth].Uvap;
      Utotal.Ugas+=U_MP[nth].Ugas;
      Utotal.Ugi+=U_MP[nth].Ugi;
      Utotal.Ugg+=U_MP[nth].Ugg;
      Utotal.Uvg+=U_MP[nth].Uvg;
      Utotal.Uvi+=U_MP[nth].Uvi;
      Utotal.Uvv+=U_MP[nth].Uvv;
    }
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
	std::vector<Atom> atomGas(int gastype);
	std::vector<vector<vector<double>>> pair_coeff;

	std::vector<Bond> bonds_v;
	std::vector<Angle> angles_v;
	std::vector<Dihedral> dihedrals_v;
	std::vector<Atom_type> atypes_v;
	std::vector<Bond_type> btypes_v;
	std::vector<Angle_type> ctypes_v;
	std::vector<Dihedral_type> dtypes_v;
	std::vector<vector<vector<double>>> pair_coeff_v;

	std::vector<vector<vector<double>>> pair_coeff_vi;
	std::vector<vector<double>> pair_coeff_vg;
	std::vector<vector<double>> pair_coeff_gi;

	Atom_type atypes_g;
	double pair_coeff_g[2];

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
