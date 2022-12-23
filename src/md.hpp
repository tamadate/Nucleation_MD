#pragma once
#include "variables.hpp"
#include "output/observer.hpp"
#include "potential/potential.hpp"
#include "PhysicalProp.hpp"
#include "flags.hpp"
#include "boundary/MBdist.hpp"

//------------------------------------------------------------------------

class MD {
	private:


  public:

	double startTime;
	int Nth;
	int calculation_number;

	int gastype;	/*1:He, 2:Ar, 3:N2*/
	int vaportype;	/*1:MeOH, 2:H2O, 3:EtOH*/
	long int step_relax;
	long int step_repre;
	long int Noftimestep;
	double p;
	double T;
	int Nof_around_gas;
	int Nof_around_vapor;
	int OBSERVE;
	double Ecoeff[3];

	double dt;
	double CUTOFF;
	double MARGIN;
	double ML2;
	double CL2;
	double d_size;
	double V;
	void setCondition(char* condfile);
	void readCondFile(char* condfile);

	long int itime;
	std::vector<long int> collisionFlagGas;
	std::vector<long int> collisionFlagVapor;
	std::vector<Potential*> InterInter;
	std::vector<Potential*> IntraInter;
	void setPotential(FLAG *flags,int mode);

	Variables *vars;
	Observer *obs;
	Physical *pp;
	FLAG *flags;
	MBdist *mbdist;
	MBdist *mbdistV;

//	General functions
	void updateInCenters(void);

//	velocity verlet
	void run_diff(char** argv);
	void verlet(void);
	void update_position(void);
	void velocity_calculation(void);
	void forceCombine(void);

//	pair list
	void updateInOut(void);
	void make_pair(void);
	void makePairXX(int type1, int type2, int cut);
	void make_pair_gasgas(void);
	void make_pairLJ(void);
	void make_pair_gasvapor(void);
	void make_pair_vaporion(void);
	void make_pairsLJCoulHybrid();
	void make_pairLJHybrid(void);
	void make_pairLJCoul(void);
	void check_pairlist(void);
  void makeDiatomicProp_in(Molecule &gasOut);
	void makePolyatomicProp_in(Molecule &vapOut);
	double margin_length;

//	initialization
	void initialization_gas(void);
  void initialization_vapor(void);
	void takeOver(void);

//	periodic
	void boundary_scaling_gas_move(void);
	void boundary_scaling_ion_move(void);
	void boundary_scaling_vapor_move(void);
	double weightFunc(double dr2);
	int loop, loop_update;	/*	current fixing time(loop) and update fixing time(loop) of out_gas for multi-timestep	*/
	double pre_ion[3];

//	analysis (calculating position and velocity of center of mass)
	void analysis_ion(void);	/*	calculation of center of ion1 and ion2, also collision judgement of collision and not collision	*/
	double gyration;

//	export
	void exportDumpOut(void);
	void exportDumpIn(void);

// related vapor sticking position
	void positionLog(void);
	int positionLogStep;
	std::vector<int> stickPositionList;
	int totalVaporIn;

/*other*/
	void output(void);
	void output_temp(double gastemp, double iontemp);
	void output_initial(void);
	void output_gas_collision(long int initime);
	void output_vapor_collision(long int initime);
	void Ovin(int i);
	void Ovout(int i);
	void display(int output_ONOFF);
	char filepath[100];
	char atomFile[100];
	char vaporFile[100];
	char takeOverFile[100];
	char vaporStickFile[100];
	char filepath_gyration[100];
	void gyration_initial(void);
	void gyration_out(MD *md2);
	string gyration_path;
	double crsq;

	void velocity_scaling(void);
	void nosehoover_ion(void);
	void nosehoover_zeta(void);
	void nosehoover_gas(void);
	void nosehoover_zeta_gas(void);
	double zeta;
	void setNVE(void);
	void setNVTion(double temp);


	double del2,CD2,rmin2;
	double totalPotential;


	MD(char* condfile,int calcNumber);
	~MD(void) {
		delete vars;
		delete obs;
		delete pp;
		delete flags;
		delete mbdist;
		delete mbdistV;
	}
	void run(char** argv);
};




//------------------------------------------------------------------------
