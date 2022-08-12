//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------
void 
MD::velocity_scaling(void) {
	double Tnow=obs->ion_temperature(vars);
	double Tp;
	Tp=300;
//	Tp=2500-2.2e-4*vars->time; 
	for (auto &a : vars->ions){
		double ratio=sqrt(Tp/Tnow);
		a.px*=ratio;
		a.py*=ratio;
		a.pz*=ratio;
	}
}

void
MD::nosehoover_zeta(void){
	double Tnow=obs->ion_temperature(vars);
	int g=vars->ions.size()*3;
	double Q_inv = 0.0001;
	vars->zeta_ion += (Tnow - pp->Tnh_ion)*g*kb_real*Q_inv*dt;
}


void
MD::nosehoover_ion(void){
	double Coeff=exp(-vars->zeta_ion*0.5*dt);
	for (auto &a : vars->ions) {
		a.px *= Coeff;
	    a.py *= Coeff;
	    a.pz *= Coeff;
	}
}

void
MD::nosehoover_zeta_gas(void){
	double Tnow=obs->gas_temperature(vars);
	int g=Nof_around_gas*3;
	double Q_inv = 0.001;
	vars->zeta_gas += (Tnow - pp->Tnh_gas)*g*kb_real*Q_inv*dt;
}

void
MD::nosehoover_gas(void){
	double Coeff=exp(-vars->zeta_gas*0.5*dt);
	for (auto &a : vars->gases) {
		a.px *= Coeff;
	    a.py *= Coeff;
	    a.pz *= Coeff;
	}
}

void
MD::setNVE(void){
	flags->velocity_scaling=0;
	flags->nose_hoover_ion=0;
	flags->nose_hoover_gas=0;
}



void 
MD::setNVTion(double temp){
	flags->velocity_scaling=0;
	flags->nose_hoover_ion=0;
	flags->nose_hoover_gas=0;
	pp->Tnh_ion=temp;
}


