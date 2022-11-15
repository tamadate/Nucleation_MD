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

struct data {
  std::vector<double> MD;
};

std::vector<double> MSDx;
std::vector<double> MSDy;
std::vector<double> MSDz;
std::vector<double> MSD;
std::vector<double> VAF;
std::vector<int> error;

double integral(std::vector<double> vaf, int blockSize){
  double sim2=0.0;
  int imax=blockSize*0.5-1;
  for(int i=1;i<imax;i++) sim2+=vaf[2*i];
  double sim3=0.0;
  imax+=1;
  for(int i=1;i<imax;i++) sim3+=vaf[2*i-1];
  return vaf[0]+vaf[blockSize-1]+2*sim2+4*sim3;
}

double slope(std::vector<double> msd, int start, int end, double dt){
  double Time=(end-start)*dt;
  return (msd[end]-msd[start])/Time;
}
