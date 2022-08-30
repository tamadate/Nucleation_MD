#pragma once
#include "struct.hpp"
#include "PhysicalProp.hpp"

//------------------------------------------------------------------------
class Variables {
public:

	Variables(void) {
		time = 0.0;
    #pragma omp parallel
    {
      #pragma omp single
      {
        Nth=omp_get_num_threads();
      }
    }
    for (int nth=0;nth<Nth;nth++){
      UionMP.push_back(0);
      UgasMP.push_back(0);
      UvapMP.push_back(0);
      UgiMP.push_back(0);
      UggMP.push_back(0);
      UvgMP.push_back(0);
      UviMP.push_back(0);
      UvvMP.push_back(0);
    }
	}

	/*variables*/
  int Nth;
	std::vector<Atom> ions;
	std::vector<Molecule> gases;
  std::vector<Molecule> vapors;
	double time;
	double zeta_ion;
	double zeta_gas;

	double Uion;
	double Ugas;
	double Uvap;
	double Ugi;
	double Ugg;
	double Uvg;
	double Uvi;
	double Uvv;

  std::vector<double> UionMP;
  std::vector<double> UgasMP;
  std::vector<double> UvapMP;
  std::vector<double> UgiMP;
  std::vector<double> UggMP;
  std::vector<double> UvgMP;
  std::vector<double> UviMP;
  std::vector<double> UvvMP;

  double tion;
  double tgas;
  double tgi;
  double tvap;
  double tvv;
  double tvg;
  double tvi;
  double tpair;


	void Uzero(void)	{
    Uion=Ugas=Uvap=Ugi=Ugg=Uvg=Uvi=Uvv=0;
    for (int nth=0;nth<Nth;nth++){
      UionMP[nth]=0;
      UgasMP[nth]=0;
      UvapMP[nth]=0;
      UgiMP[nth]=0;
      UggMP[nth]=0;
      UvgMP[nth]=0;
      UviMP[nth]=0;
      UvvMP[nth]=0;
    }
  }
  void Ucombine(void)	{
    for (int nth=0;nth<Nth;nth++){
      Uion+=UionMP[nth];
      Uvap+=UvapMP[nth];
      Ugas+=UgasMP[nth];
      Ugi+=UgiMP[nth];
      Ugg+=UggMP[nth];
      Uvg+=UvgMP[nth];
      Uvi+=UviMP[nth];
      Uvv+=UvvMP[nth];
    }
  }
  void tzero(void)	{tion=tgas=tvap=tgi=tvv=tvg=tvi=tpair=0;}
	double Usum(void)	{return Uion+Ugas+Uvap+Ugi+Ugg+Uvg+Uvi+Uvv;}

	std::vector<int> gas_in;	/*	gas list around ion1	*/
	std::vector<int> gas_out;	/*	gas list far from ion1	*/
	std::vector<int> vapor_in;	/*	gas list around ion1	*/
	std::vector<int> vapor_out;	/*	gas list far from ion1	*/

	/*vectors for potential calculation*/
  std::vector<Pair> ion_pairs;
	std::vector<Pair> pairs_gi;	/*	gas-ion interaction pair list	*/
	std::vector<Pair> pairs_gv;	/*	gas-vapor interaction pair list	*/
	std::vector<Pair> pairs_gg;	/*	gas-gas interaction pair list	*/
	std::vector<Bond> bonds;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;
	std::vector<Atom_type> atypes;
	std::vector<Bond_type> btypes;
	std::vector<Angle_type> ctypes;
	std::vector<Dihedral_type> dtypes;
	std::vector<int> molecules;
	std::vector<Atom> atomVapor;
	std::vector<Bond> bondMeOH(void);
	std::vector<Angle> angleMeOH(void);
	std::vector<Dihedral> dihedralMeOH(void);
	std::vector<Atom> atomGas(void);
	std::vector<vector<vector<double>>> pair_coeff;
	double bornCoeff[2][2][5];

	/*initialization and export to dump file*/
	void read_initial(char* infile);
	void set_initial_velocity(Physical *pp);
	double totalPotential;
	double totalVirial;

	void ROTATION(double &X, double &Y, double &Z, double A, double B, double C, double x, double y, double z){
		X = cos(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+sin(C)*(y*cos(A)+z*sin(A));
		Y = -sin(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+cos(C)*(y*cos(A)+z*sin(A));
		Z = x*sin(B)-y*sin(A)*cos(B)+z*cos(A)*cos(B);
	}

	std::vector<Atom> makeAtomMeOH(void);
	std::vector<Atom> makeAtomTIP3P(void);


private:
};
//------------------------------------------------------------------------
