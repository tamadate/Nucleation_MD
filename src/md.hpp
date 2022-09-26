#pragma once
#include "variables.hpp"
#include "observer.hpp"
#include "potential.hpp"
#include "PhysicalProp.hpp"
#include "flags.hpp"
#include "MBdist.hpp"
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

//	vectors for pairlist
	double margin_length;

//	General functions
	void updateGasinCenters(void);
	void updateVaporinCenters(void);

//	velocity verlet
	void run_diff(char** argv);
	void verlet(void);
	void update_position(void);
	void velocity_calculation(void);
  void update_position_constrained(void);
  void update_velocity_constrained(void);

//	pair list
	void update_vapor_in(void);
	void update_gas_in(void);
	void make_pair(void);
	void make_pair_gasgas(void);
	void make_pair_gasion(void);
	void make_pair_gasvapor(void);
	void make_pair_vaporion(void);
	void make_pair_vaporvapor();
	void check_pairlist(void);
  void makeDiatomicProp_in(int i);
  void makeDiatomicProp_out(int i);
	void makePolyatomicProp_in(int i);
	void makePolyatomicProp_out(int i);

//	initialization
	void initialization_gas(void);
  void initialization_vapor(void);

//	periodic
	void periodic(void);	/*	periodic condition for gas_in	*/
	void boundary_scaling_gas_move(void);
	void boundary_scaling_ion_move(void);
	void boundary_scaling_vapor_move(void);
	int loop, loop_update;	/*	current fixing time(loop) and update fixing time(loop) of out_gas for multi-timestep	*/
	double pre_ion[3];

//	analysis (calculating position and velocity of center of mass)
	void analysis_gas(void);	/*	calculation of center of ion1 and ion2, also collision judgement of collision and not collision	*/
	void analysis_ion(void);	/*	calculation of center of ion1 and ion2, also collision judgement of collision and not collision	*/
	double ion_r[3];
	double ion_v[3];	/*	center of ion1 and ion2	*/
	double ion_f[3];
	double gas_r[3];
	double gas_v[3];	/*	center of ion1 and ion2	*/
	double gyration;

//	export
	void export_dump(void);
	void export_dump_close(void);

// related vapor sticking position
	void positionLog(void);
	int positionLogStep;
	std::vector<int> stickPositionList;

/*other*/
	void output(void);
	void output_gas(void);
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
	char vaporStickFile[100];
	char filepath_gyration[100];
	void fix_cell_center(void);
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
	~MD(void);
	void run(char** argv);
	int yesno;	/*	flag for collision or not collision	*/
};




//------------------------------------------------------------------------
