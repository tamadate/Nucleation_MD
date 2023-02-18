#pragma once
#include "../variables.hpp"
//------------------------------------------------------------------------
class Observer {
public:
	Variables *vars;
	double startTime;
	double V;
	double d_size;
	const double coeff=kb_real_inv/1.5;
	double Kin_g;
	double Kout_g;
	double K_g;
	double Tin_g;
	double Tout_g;
	double T_g;

	double Kin_v;
	double Kout_v;
	double K_v;
	double Tin_v;
	double Tout_v;
	double T_v;

	double Kion;
	double Tion;

	int Ngas;
	int Nion;
	int Nvap;

	void computeGasProps(void);
	void computeIonProps(void);
	void computeVaporProps(void);
	double potential_energy() {return vars->totalPotential;}
	double pressure(std::vector<Pair> &pairs, double Treal, double virial,double p,double T);
	double total_energy(std::vector<Pair> &pairs) {return K_g + potential_energy();}

	void outputDump(void);
	void outputDumpClose(void);
	void outputIonCenter(void);
	void outputGasCenter(void);
	void Ovin(int i, double time);
	void Ovout(int i, double time);
	void display(void);

	char fileVaporIn[100];
	char fileVaporOut[100];
	char fileKinetic[100];
	char filePotential[100];
	char fileIonCenter[100];
	char fileGasCenter[100];
	char *fileDump;

	bool dump_fix;
	Observer(Variables *VARS, int calculation_number);
};
//------------------------------------------------------------------------
