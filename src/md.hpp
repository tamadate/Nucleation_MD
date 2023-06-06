#pragma once
#include "variables.hpp"
#include "output/observer.hpp"
#include "potential/potential.hpp"
#include "PhysicalProp.hpp"
#include "pairlist/MBdist.hpp"
#include "MDcondition.hpp"
#include "thermostat/thermostat.hpp"
#include "vaporUptake/collision_inout/collision.hpp"
#include "sampling/gyration/gyration.hpp"
#include "sampling/ionCOM/ionCenter.hpp"
#include "vaporUptake/stickPosition/stickPosition.hpp"
#include "vaporUptake/collision_trajectory/trajectory.hpp"
#include "functions.hpp"

class MD {
private:


public:
	int gastype;	// gas type
	long int itime; // number of iteration
	
	std::vector<Potential*> InterInter; // inter-molecular interactions
	std::vector<Potential*> IntraInter;	// intra-molecular interactiops

	Variables *vars;	// variables class
	Observer *obs;		// observer class
	Physical *pp;		// physical property class
	MDcondition *con;	// calculation condition class
	MBdist *mbdist;		// Maxwell-Boltzumann distribution class for gas moleucle
	MBdist *mbdistV;	// Maxwell-Boltzumann distribution class for vapor molecule
	std::vector<Functions*> funcs; // funciton class array

// start MD simulation
	void run(char** argv);	

//	velocity verlet
	void verlet(void);					// main function
	void update_position(void); 		// update position in Verlet
	void velocity_calculation(void);	// update velocity in Verlet

//	pair list
	void update_vapor_in(void);			// update list of all-atom vapor (vars->vapor_in)
	void update_gas_in(void);			// update list of all-atom gas (vars->gas_in)
	void make_pair(void);				// update pair list
	void check_pairlist(void);			// check if pair list is updated or not
	int loop;							// number of iterations from last pair-list updating
	int loop_update;					// number of iterations to update pair-list
	/* 
		When the molecules are transfered either from mono-atomic to poly-atomic molecule or the other direction,
		it is necessary to modify some matrixes in the code (you can find detail in below four functions).
		makeDiatomicProp is for gas molecule and makePolyatomicProp is for vapor molecule.
		Signs "in" and "out" indicate the case entering to and exiting from the all-atom region, respectively.
	*/
	void makeDiatomicProp_in(int i);	
	void makeDiatomicProp_out(int i);
	void makePolyatomicProp_in(int i);
	void makePolyatomicProp_out(int i);
	void updateInCenters(void);

//	initialization
	void initialization_gas(void);		// initialization of gas molecules
	void initializatIonVapor(void);		// initialization of ion and vapor molecules
	void readCondFile(char* condfile); // read conditions file

//	periodic
	void periodic(void);	//	periodic condition for gas_in	
	/* 
		Below 3 functions are velocity rescaling periodic bnoundary conditions.
		It resample the gas and vapor molecule velocites when it re-enter to the simulation domain 
		through periodic B.C.
		There are two types which sample the velocity from normal MB and modified MB distributions
	*/
	void boundary_scaling_gas_move(void);	// modified MB-distribution for gas molecule
	void boundary_scaling_vapor_move(void);	// modified MB-distribution for vapor molecule
	void boundary_scaling_ion_move(void);	// modified MB-distribution for gas molecule
	double pre_ion[3];

//	analysis (calculating position and velocity of center of mass)
	void getGasCenterProp(void);	// For gas molecules
	void getIonCenterProp(void);	// For ion
	// distance of atoms in the ion from ion COM 
	double distFromIonCenter(Molecule &mol, double &dx, double &dy, double &dz); 	
	// distance of atoms in the ion from "previous" ion COM 
	double distFromIonPreCenter(Molecule &mol, double &dx, double &dy, double &dz);


/*other*/
	char filepath[100]; // general file path
	char atomFile[100]; // file path for ion
	char vaporFile[100]; // file path for vapor


	MD(char* condfile,int calcNumber);
	~MD(void);
};




//------------------------------------------------------------------------
