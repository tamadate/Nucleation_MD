#pragma once

#define _USE_MATH_DEFINES
#include <algorithm>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <random>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "struct.hpp"
#include "flags.hpp"
//#pragma GCC target("avx2")

using namespace std;
//------------------------------------------------------------------------

const double Nw = 6.02e23;
const double kb = 1.38e-23;     /* boltzmann constant */
const double kb_real = kb*Nw/4184.0;     /* boltzmann constant with real unit kcal/molK*/
const double e0 = 8.85e-12;     /* permittivity of vacuum */
const double e = 1.602e-19;     /* elementary charge */

const double MN2=28.0, MHe=4.0026, MAr=39.948;
const double MMeOH=32.0, MH2O=18.0, MEtOH=46.0;
const double mN2=MN2/Nw/1000.0;
const double myuHe=1.47e-5;
const double alphaHe=0.208, alphaN2=1.7*0.5, alphaAr=1.664;
const double SAr=114, TrefAr=273.0, myuAr=2.125e-5; // REF
const double SN2=107, TrefN2=273.0, myuN2=1.663e-5; // REF

//------------------------------------------------------------------------

/*******************Coeff.************************/
const double Rinter=100;	   // boundary of AA and CG models
const double BoundL=10;     // Thickness of overlapping region
const double BoundMergin=10;   // Mergin of overlapping region
const double RAA=Rinter+BoundL;   // Mergin of overlapping region
const double RAA2=RAA*RAA;
const double RCG=Rinter-BoundL;   // Mergin of overlapping region
const double RCG2=RCG*RCG;
const double RI2=(RAA+BoundMergin)*(RAA+BoundMergin); // AA outside boundary radius
const double RO2=(RCG-BoundMergin)*(RCG-BoundMergin); // CG hollow boundary radius
const double cal=4187;
const double qqrd2e=e*e/4.0/M_PI/e0*Nw*1e10/4184.0;
const double Cpress=1e30/Nw*4184;
const double eV_to_kcalmol=23.061;
const double real_to_kcalmol=10000/4.184;	/*	A^2 g fs^-2 mol^-1 to kcal/mol*/
const double kb_real_inv = 1/kb_real;     /* boltzmann constant with real unit kcal/molK*/
const int AA=00000001;
const int CG=00000010;
const int AACG=00000011;

void adjust_periodic(double &dx, double &dy, double &dz, double d_size);
