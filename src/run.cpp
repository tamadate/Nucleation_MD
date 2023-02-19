//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*
	- Diffusion coefficient calculation
	- MD simulation performed for (step_relax) steps as a thermal
	relaxation, and main calculation run for (Noftimestep) steps.
	- Among of two calculations, a calculation is re-initialized.
	- Properties (pressure, temperature, energy etc...) are observed
	each (OBSERVE) steps.
*/
/////////////////////////////////////////////////////////////////////
void
MD::run(char** argv) {
	/****Thermal relaxation****/
	for (itime=0; itime < con->step_relax; itime++) {
		vars->time=itime*dt;
		if (itime%obs->OBSERVE==0) {
			obs->display();
			obs->outputDumpClose();
		}
    	verlet();
	}

	/*
	****Reinitialization****
	- Reset time (time -> 0).
	- Re-initializating the ion and gas positions and velocities.
	Position -> 0 elta, Traslational velocity -> 0
	- Set ion's center of mass (maybe -> 0), make pair list for initial
	step of simulation, reset the margine size.
	*/
	for (auto &a : vars->ions) a.qx-=vars->IonX[0], a.qy-=vars->IonX[1], a.qz-=vars->IonX[2];
	for (auto &a : vars->gases) a.qx-=vars->IonX[0], a.qy-=vars->IonX[1], a.qz-=vars->IonX[2];
	getIonCenterProp();
	make_pair();
	delete thermo;
	thermo = new Thermostat();

	/****Main simulaiton****/
	for (itime=0; itime<con->Noftimestep; itime++) {
		vars->time=itime*dt;
		if (itime%obs->OBSERVE==0) {
			flags->eflag=true;
			vars->Uzero();
		}
		if (con->positionLogStep>0){
			if(itime%con->positionLogStep==0) positionLog();
		}
		if(itime%con->logger==0){
			getGasCenterProp();
			obs->outputIonCenter();
			//obs->outputGasCenter();
		}
		if (itime%obs->OBSERVE==0) {
			obs->display();
			obs->outputDumpClose();
			flags->eflag=false;
		}
		verlet();
	}
}
