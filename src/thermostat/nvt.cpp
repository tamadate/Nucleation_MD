//------------------------------------------------------------------------
#include "thermostat.hpp"
//------------------------------------------------------------------------
void
ThermostatVscale::tempControl(double dt) {
	obs->computeIonProps();
	for (auto &a : vars->ions){
		double ratio=sqrt(Ttarget/obs->Tion);
		a.px*=ratio;
		a.py*=ratio;
		a.pz*=ratio;
	}
}

void
ThermostatNH::computeZeta(double dt){
	obs->computeIonProps();
	int g=vars->ions.size()*3;
	zeta += (obs->Tion - Ttarget)*g*kb_real*Q_inv*dt;
}


void
ThermostatNH::tempControl(double dt){
	double Coeff=exp(-zeta*0.5*dt);
	for (auto &a : vars->ions) {
		a.px *= Coeff;
	    a.py *= Coeff;
	    a.pz *= Coeff;
	}
}

