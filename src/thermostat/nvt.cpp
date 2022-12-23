//------------------------------------------------------------------------
#include "../md.hpp"
//------------------------------------------------------------------------
void
MD::velocity_scaling(void) {
	obs->computeProps(vars,0);
	double Tp;
	Tp=300;
//	Tp=2500-2.2e-4*vars->time;
	for (auto &a : vars->Molecules[0].inAtoms){
		double ratio=sqrt(Tp/obs->Tin[0]);
		a.px*=ratio;
		a.py*=ratio;
		a.pz*=ratio;
	}
}


void
MD::nosehoover_zeta(void){
	obs->computeProps(vars,0);
	int g=vars->Molecules[0].inAtoms.size()*3;
	double Q_inv = 0.0001;
	vars->zeta_ion += (obs->Tin[0] - pp->Tnh_ion)*g*kb_real*Q_inv*dt;
}


void
MD::nosehoover_ion(void){
	double Coeff=exp(-vars->zeta_ion*0.5*dt);
	for (auto &a : vars->Molecules[0].inAtoms) {
		a.px *= Coeff;
    a.py *= Coeff;
    a.pz *= Coeff;
	}
}

void
MD::nosehoover_zeta_gas(void){
	obs->computeProps(vars,1);
	int g=Nof_around_gas*3;
	double Q_inv = 0.001;
	vars->zeta_gas += (obs->Tout[1] - pp->Tnh_gas)*g*kb_real*Q_inv*dt;
}

void
MD::nosehoover_gas(void){
	double Coeff=exp(-vars->zeta_gas*0.5*dt);
	for (auto i : vars->MolID[1]) {
		vars->Molecules[i].px *= Coeff;
    vars->Molecules[i].py *= Coeff;
    vars->Molecules[i].pz *= Coeff;
	}
}


void
MD::setNVTion(double temp){
	flags->velocity_scaling=0;
	flags->nose_hoover_ion=0;
	flags->nose_hoover_gas=0;
	pp->Tnh_ion=temp;
}
