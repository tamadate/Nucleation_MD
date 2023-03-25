#pragma once
#include "../functions.hpp"

class NVT_vscale : public Functions {
	private:
		string statName="NVT velocity scaling";
	public:
		double Ttarget;
		virtual void printName(void) {cout<<statName<<endl;}
		void postLoop(void){
			obs->computeIonProps();
			for (auto &a : vars->ions){
				double ratio=sqrt(Ttarget/obs->Tion);
				a.px*=ratio;
				a.py*=ratio;
				a.pz*=ratio;
			}
		};
		void initial(){};
		NVT_vscale(double T, int Step){
			Ttarget=T;
			type=1;
			step=Step;
		};
		~NVT_vscale(){};
};
