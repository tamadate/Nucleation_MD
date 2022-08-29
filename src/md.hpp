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

	long int itime;
	std::vector<long int> collisionFlagGas;
	std::vector<long int> collisionFlagVapor;
	std::vector<Potential*> InterInter;
	std::vector<Potential*> IntraInter;
	void setPotential(FLAG *flags);

	Variables *vars;
	Observer *obs;
	Physical *pp;
	FLAG *flags;
	MBdist *mbdist;
	MBdist *mbdistV;

//	vectors for pairlist
	double margin_length;

//	velocity verlet
	void run_diff(char** argv);
	void verlet(void);
	void update_position(void);
	void velocity_calculation(void);
  void update_position_constrained(void);
  void update_velocity_constrained(void);
	void forceCombine(void);

	//	pair list
	void update_vapor_in(void);
	void update_gas_in(void);
	void make_pair(void);
	void make_pair_gasgas(void);
	void make_pair_gasion(void);
	void make_pair_gasvapor(void);
	void make_pair_vaporion(void);
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



		MD(char* condfile);
		~MD(void);
		void run(char** argv);
		int yesno;	/*	flag for collision or not collision	*/
};




//------------------------------------------------------------------------
