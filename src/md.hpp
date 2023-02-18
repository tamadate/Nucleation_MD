#pragma once
#include "variables.hpp"
#include "output/observer.hpp"
#include "potential/potential.hpp"
#include "PhysicalProp.hpp"
#include "flags.hpp"
#include "pairlist/MBdist.hpp"
#include "MDcondition.hpp"
#include "thermostat/thermostat.hpp"
//------------------------------------------------------------------------

class MD {
	private:


	public:

	int Nth;
	int calculation_number;

	int gastype;	/*1:He, 2:Ar, 3:N2*/
	
	double dt;
	long int itime;
	
	std::vector<long int> collisionFlagGas;
	std::vector<long int> collisionFlagVapor;
	std::vector<Potential*> InterInter;
	std::vector<Potential*> IntraInter;

	Variables *vars;
	Observer *obs;
	Physical *pp;
	FLAG *flags;
	MDcondition *con;
	Thermostat *thermo;
	MBdist *mbdist;
	MBdist *mbdistV;

//	General functions

//	velocity verlet
	void run(char** argv);
	void verlet(void);
	void update_position(void);
	void velocity_calculation(void);

//	pair list
	void update_vapor_in(void);
	void update_gas_in(void);
	void make_pair(void);
	void check_pairlist(void);
	void makeDiatomicProp_in(int i);
	void makeDiatomicProp_out(int i);
	void makePolyatomicProp_in(int i);
	void makePolyatomicProp_out(int i);
	void updateInCenters(void);

//	initialization
	void initialization_gas(void);
	void initializatIonVapor(void);
	void readCondFile(char* condfile);

//	periodic
	void periodic(void);	/*	periodic condition for gas_in	*/
	void boundary_scaling_gas_move(void);
	void boundary_scaling_ion_move(void);
	void boundary_scaling_vapor_move(void);
	int loop, loop_update;	/*	current fixing time(loop) and update fixing time(loop) of out_gas for multi-timestep	*/
	double pre_ion[3];

//	analysis (calculating position and velocity of center of mass)
	void getGasCenterProp(void);	/*	calculation of center of ion1 and ion2, also collision judgement of collision and not collision	*/
	void getIonCenterProp(void);	/*	calculation of center of ion1 and ion2, also collision judgement of collision and not collision	*/
	double distFromIonCenter(Molecule &mol, double &dx, double &dy, double &dz);
	double distFromIonPreCenter(Molecule &mol, double &dx, double &dy, double &dz);
	double gyration;


// related vapor sticking position
	void positionLog(void);
	std::vector<int> stickPositionList;

/*other*/

	char filepath[100];
	char atomFile[100];
	char vaporFile[100];
	char vaporStickFile[100];
	char filepath_gyration[100];
	void gyration_initial(void);
	void gyration_out(MD *md2);
	string gyration_path;
	double crsq;

	double totalPotential;

	MD(char* condfile,int calcNumber);
	~MD(void);
};




//------------------------------------------------------------------------
