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

using namespace std;


const double Nw = 6.02e23;      // Avogadro number
const double kb = 1.38e-23;     // Boltzmann constant 
const double kb_real = kb*Nw/4184.0;     // Boltzmann constant with real unit kcal/molK
const double e0 = 8.85e-12;     // permittivity of vacuum 
const double e = 1.602e-19;     // elementary charge 

// gas molecule aconstant (He, Ar, and N2)
const double MN2=28.0, MHe=4.0026, MAr=39.948;  // molecular weight
const double alphaHe=0.208, alphaN2=1.7*0.5, alphaAr=1.664; // polarizability
const double SAr=114, TrefAr=273.0, myuAr=2.125e-5; // constant for Sutherland's equation
const double SN2=107, TrefN2=273.0, myuN2=1.663e-5; // constant for Sutherland's equation

// coefficients and conversiont factors
const double Rinter=200;	            // Radius of all-atom region, R
const double RI2=Rinter*Rinter;         // square of R
const double cal=4187;                  // conversion factor between kcal and J
const double qqrd2e=e*e/4.0/M_PI/e0*Nw*1e10/4184.0; // coefficient of electrical force/potential calculation
const double Cpress=1e30/Nw*4184;       // conversion factor for pressure unit conversion
const double eV_to_kcalmol=23.061;      // conversion factor between eV and kcal/mol
const double real_to_kcalmol=10000/4.184;	//	conversion factor from A^2 g fs^-2 mol^-1 to kcal/mol
const double kb_real_inv = 1/kb_real;     // Boltzmann constant with real unit kcal/molK

void adjust_periodic(double &dx, double &dy, double &dz, double d_size);
