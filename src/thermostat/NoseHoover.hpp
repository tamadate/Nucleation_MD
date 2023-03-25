#pragma once
#include "../functions.hpp"

class NVT_nh : public Functions {
	private:
		string statName="NVT Nose-Hoover";
	public:
		double Ttarget;
		double zeta;
		double Q_inv;
		virtual void printName(void) {cout<<statName<<endl;}
		void postLoop(void){
			double Coeff=exp(-zeta*0.5*con->dt);
			for (auto &a : vars->ions) {
				a.px *= Coeff;
				a.py *= Coeff;
				a.pz *= Coeff;
			}
		};
		void postPosition(void){
			obs->computeIonProps();
			int g=vars->ions.size()*3;
			zeta += (obs->Tion - Ttarget)*g*kb_real*Q_inv*con->dt;
		};
		void initial(){};
		NVT_nh(double T, int Step){
			Ttarget=T;
			type=1;
			Q_inv=0.0001;
			step=Step;
		};
		~NVT_nh(){};
};
