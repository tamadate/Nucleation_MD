#include "../md.hpp"

void
Observer::computeGasProps(void){
	Kin_g=0;
	double Nin=0;
	Kout_g=0;
	double Nout=0;

	Molecule *gases = vars->gases.data();
	for(auto &I : vars->gas_in){
		for (auto &ag : gases[I].inAtoms){
			Kin_g += ag.px * ag.px * ag.mass;
			Kin_g += ag.py * ag.py * ag.mass;
			Kin_g += ag.pz * ag.pz * ag.mass;
			Nin++;
		}
	}
	Kin_g*= (0.5 * real_to_kcalmol);
	Tin_g=Kin_g/double(Nin)*coeff;

	for(auto i : vars->gas_out){
		Kout_g += gases[i].px * gases[i].px * gases[i].mass;
		Kout_g += gases[i].py * gases[i].py * gases[i].mass;
		Kout_g += gases[i].pz * gases[i].pz * gases[i].mass;
		Nout++;

	}
	Kout_g*= (0.5 * real_to_kcalmol);
	Tout_g=Kout_g/double(Nout)*coeff;

	K_g=Kin_g+Kout_g;
	Ngas=Nin+Nout;
	T_g=K_g/double(Ngas)*coeff;
};

void
Observer::computeIonProps(void){
	Kion=0;
	Nion=0;
	for (auto &a : vars->ions) {
		Kion += a.px * a.px * a.mass;
		Kion += a.py * a.py * a.mass;
		Kion += a.pz * a.pz * a.mass;
		Nion++;
	}
	Kion*= (0.5 * real_to_kcalmol);
	Tion=Kion/double(Nion)*coeff;
};

void
Observer::computeVaporProps(void){
	Kin_v=0;
	double Nin=0;
	Kout_v=0;
	double Nout=0;

	Molecule *vapors = vars->vapors.data();
	vars->vapor_in.size();
	for(auto &I : vars->vapor_in){
		for (auto &av : vapors[I].inAtoms){
			Kin_v += av.px * av.px * av.mass;
			Kin_v += av.py * av.py * av.mass;
			Kin_v += av.pz * av.pz * av.mass;
			Nin++;
		}
	}
	Kin_v*= (0.5 * real_to_kcalmol);
	if(Nin>0)	Tin_v=Kin_v/double(Nin)*coeff;
	else Tin_v=0;

	for(auto i : vars->vapor_out){
		Kout_v += vapors[i].px * vapors[i].px * vapors[i].mass;
		Kout_v += vapors[i].py * vapors[i].py * vapors[i].mass;
		Kout_v += vapors[i].pz * vapors[i].pz * vapors[i].mass;
		Nout++;
	}
	Kout_v*= (0.5 * real_to_kcalmol);
	Tout_v=Kout_v/double(Nout)*coeff;

	K_v=Kin_v+Kout_v;
	Nvap=Nin+Nout;
	T_v=K_v/double(Nvap)*coeff;
};