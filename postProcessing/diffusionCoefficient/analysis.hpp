#pragma once
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>
#include <stdlib.h>
using namespace std;

const double e=1.602e-19;
const double kb=1.38e-23;
const double T=300.0;
const double P=1e5;
const double Nw=6.02e23;
const double mi=1046.2/1000.0/Nw; // mass of ion
const double mg=28/1000.0/Nw;     // mass of gas
const double mred=1/(1/mg+1/mi);  // reduced mass (ion-gas)
const int sample_step=1000;
const int start=500;
const int dstep=1;
const int output=10000;           // output interval of simulation
const int total_step=0.5e9/output;// total 
const int skip=0e9/output;
const int N=(total_step-sample_step)/dstep;
const double dt=output*1e-15;
const int pl=1;

double t[sample_step], MSD[sample_step], MSDx[sample_step], MSDy[sample_step], MSDz[sample_step], VAF[sample_step];
double velocityAve[3][sample_step];
double data[total_step+1][8], data_gas[total_step+1][8];
