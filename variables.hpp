#pragma once
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include "CalculationCond.hpp"
#include "PhysicalProp.hpp"

using namespace std;
//------------------------------------------------------------------------
struct Atom_type {
 	int type;
	double mass;
	double coeff1, coeff2;
};
//------------------------------------------------------------------------
struct Bond_type {
 	int type;
	double coeff[2];
};
//------------------------------------------------------------------------
struct Angle_type {
 	int type;
	double coeff[2];
};
//------------------------------------------------------------------------
struct Dihedral_type {
 	int type, multi;
	double coeff[15];
};
//------------------------------------------------------------------------
struct Ion {
 	double qx, qy, qz;
	double px, py, pz;
	double charge;
	int type;
	int id;
	double fx, fy, fz;
	double mass;
	int ix, iy, iz;
};
//------------------------------------------------------------------------
struct Gas {
 	double qx, qy, qz;
	double px, py, pz;
	double charge;
	int type;
	int id;
	double fx, fy, fz;
	double mass;
	int ix, iy, iz;
};
//------------------------------------------------------------------------
struct Atom {
 	double qx, qy, qz;
	double px, py, pz;
	double charge;
	int type;
	int id;
	double fx, fy, fz;
	double mass;
	int ix, iy, iz;
};

//------------------------------------------------------------------------
struct Bond {
	int atom1, atom2, type;
};
//------------------------------------------------------------------------
struct Angle {
	int atom1, atom2, atom3, type;
};
//------------------------------------------------------------------------
struct Dihedral {
	int atom1, atom2, atom3, atom4, type;
};
//------------------------------------------------------------------------
struct Pair{
	int i,j;
};
struct Pair_many{
	int i;
	std::vector<int> j;	/*	pair list	*/
};

//------------------------------------------------------------------------
struct N2_relative{
	double x,y,z,vx1,vy1,vz1,vx2,vy2,vz2;
};
//------------------------------------------------------------------------
struct Molecule {
 	double qx, qy, qz;
	double px, py, pz;
	double mass;
	int ix, iy, iz;
	std::vector<Atom> inAtoms;
	std::vector<Bond> bonds;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;
	int inFlag;
};
//------------------------------------------------------------------------
class Variables {
public:

	Variables(void) {
		time = 0.0;
	}

	/*variables*/
	std::vector<Ion> ions;
	std::vector<Gas> gases;
    std::vector<Molecule> vapors;
	double time;
	double zeta_ion;
	double zeta_gas;

	std::vector<int> gas_in;	/*	gas list around ion1	*/
	std::vector<int> gas_out;	/*	gas list far from ion1	*/
	std::vector<int> vapor_in;	/*	gas list around ion1	*/
	std::vector<int> vapor_out;	/*	gas list far from ion1	*/

	/*vectors for potential calculation*/
    std::vector<Pair> ion_pairs;
	std::vector<Pair> pairs_gi;	/*	gas-ion interaction pair list	*/
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
	std::vector<vector<vector<double>>> pair_coeff;
	double bornCoeff[2][2][5];

	/*add to vectors*/
	void add_ions(int id, int type, double x, double y, double z, double vx, double vy, double vz, double fx, double fy, double fz, double charge, double mass);
	void add_gases(int id, int type, double x, double y, double z, double vx, double vy, double vz, double fx, double fy, double fz, double charge, double mass);
    void add_vapors(int id, int type, double x, double y, double z, double vx, double vy, double vz, double mass);
	void add_bonds(int atom1, int atom2, int type);
	void add_angles(int atom1, int atom2, int atom3, int type);
	void add_dihedrals(int atom1, int atom2, int atom3, int atom4, int type);
	void add_atype(int type, double mass, double coeff1, double coeff2);

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
