#pragma once
#include "variables.hpp"
#include "output/observer.hpp"
#include "potential/potential.hpp"
#include "PhysicalProp.hpp"
#include "pairlist/MBdist.hpp"
#include "MDcondition.hpp"
#include "thermostat/thermostat.hpp"
#include "sampling/collision/collision.hpp"
#include "sampling/gyration/gyration.hpp"
#include "sampling/ionCenter/ionCenter.hpp"
#include "sampling/stickPosition/stickPosition.hpp"
#include "sampling/trajectory/trajectory.hpp"
#include "functions.hpp"
//------------------------------------------------------------------------

class MD {
private:


public:
	int gastype;	/*1:He, 2:Ar, 3:N2*/
	long int itime;
	
	std::vector<Potential*> InterInter;
	std::vector<Potential*> IntraInter;

	Variables *vars;
	Observer *obs;
	Physical *pp;
	MDcondition *con;
	MBdist *mbdist;
	MBdist *mbdistV;
	std::vector<Functions*> funcs;

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


/*other*/
	char filepath[100];
	char atomFile[100];
	char vaporFile[100];
	double crsq;

	double totalPotential;

	MD(char* condfile,int calcNumber);
	~MD(void);
};




//------------------------------------------------------------------------
