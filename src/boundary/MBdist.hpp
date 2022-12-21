#pragma once
#include "../constants.hpp"
#define N 5000000
#define M 1000

//------------------------------------------------------------------------
class MBdist {
	private:
	public:
		int number;
		std::vector<double> vflux;	/*	seed of vf(v)/c	*/
		int* dis;	/*	number of gas molecuole v(dv*i)<v<v(dv*(i+1)) */

		void makeWeightedMB(double meanVel, double mass, double T){
			number=0;
			double vmax=5.0*meanVel;	/*	max velocity of gas molecule in flux velocity distribution */
			double dv=vmax/M;	/*	resolution of flux velocity distribution */
			dis = new int[M];
			vector<double>().swap(vflux);
			for(int i=0; i<M; i++) dis[i]=0;
			random_device seed;
			mt19937 mt(seed());
			uniform_real_distribution<double> v(0,vmax);
			double kbT=kb*T;
			double kbT2=kbT*kbT;
			double coeffA=0.5*mass*mass/kbT2;
			double coeffB=-0.5*mass/kbT;
			for(int i=0; i<N; i++) {
				double vgen=v(mt);
				int idist=vgen/dv;
				double x=(idist+0.5)*dv;
				double x2=x*x;
				double PDFweighted=coeffA*x2*x*exp(coeffB*x2);
				if(dis[idist]<PDFweighted*2000000){
					vflux.push_back(vgen);
					dis[idist]++;
				}
			}
//			FILE*f=fopen("weightedDist.dat", "w");
//			for(int i=0; i<M; i++) {fprintf(f, "%e %d \n", (i+0.5)*dv, dis[i]);}
//			fclose(f);
			random_shuffle(vflux.begin(), vflux.end());
			delete [] dis;
		};

		MBdist(){};
		~MBdist(){};
};
