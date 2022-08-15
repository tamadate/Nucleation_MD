#pragma once

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <cstring>
#include <string>
#include <algorithm>

using namespace std;
extern int gastype;	/*1:He, 2:Ar, 3:N2*/
extern int vaportype;	/*1:MeOH, 2:H2O, 3:EtOH*/
extern long int Noftimestep;
extern double p;
extern double T;
extern int calculation_number;
extern int Nof_around_gas;
extern int Nof_around_vapor;
extern int OBSERVE;
extern double d_size;
extern double V;
extern double Ex,Ey,Ez;


//------------------------------------------------------------------------
extern double dt;
extern double CUTOFF;
extern double MARGIN;
extern double MARGIN_PEG;
extern double ML2;
extern double ML2_V;
extern double CL2;
const double sqCollsionDistanceGasIon=5*5;
const double sqCollsionDistanceVaporIon=5*5;
//------------------------------------------------------------------------

const double Nw = 6.02e23;
const double kb = 1.38e-23;     /* boltzmann constant */
const double kb_real = kb*Nw/4184.0;     /* boltzmann constant with real unit kcal/molK*/
const double e0 = 8.85e-12;     /* permittivity of vacuum */
const double e = 1.602e-19;     /* elementary charge */
const double MN2=28.0, MHe=4.0026;
const double Rinter=100;	// Radius of interaction area
const double RI2=Rinter*Rinter;


void adjust_periodic(double &dx, double &dy, double &dz);
//------------------------------------------------------------------------

const double cal=4187;

extern double myu, D0NH4, D0NO2, DNH4, DNO2, Mgas, mNH4_gas, mNO2_gas;

/*******************Coeff.************************/
const double qqrd2e=e*e/4.0/M_PI/e0*Nw*1e10/4184.0;
const double Cpress=1e30/Nw*4184;
const double eV_to_kcalmol=23.061;
const double real_to_kcalmol=10000/4.184;	/*	A^2 g fs^-2 mol^-1 to kcal/mol*/
const double kb_real_inv = 1/kb_real;     /* boltzmann constant with real unit kcal/molK*/

/*******************Function************************/
void SET_CONDITION();
