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
  	const int logger=10000;
	for (auto &a : IntraInter) {a->printName();}
	for (auto &a : InterInter) {a->printName();}
	for (itime=0; itime < step_relax; itime++) {
		if (itime%OBSERVE==0) {
			obs->display();
			obs->export_dump_close();
			vars->time=(itime*dt);
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
	vars->time=0;
	itime=0;
	for (auto &a : vars->ions) a.qx-=vars->ion_r[0], a.qy-=vars->ion_r[1], a.qz-=vars->ion_r[2];
	for (auto &a : vars->gases) a.qx-=vars->ion_r[0], a.qy-=vars->ion_r[1], a.qz-=vars->ion_r[2];
	getIonCenterProp();
	make_pair();
	margin_length = MARGIN;
	setNVE();

/****Main simulaiton****/
	for (itime=0; itime<Noftimestep; itime++) {
    if (positionLogStep>0){
		if(itime%positionLogStep==0) positionLog();
    }
    if(itime%logger==0){
		getGasCenterProp();
		obs->outputIonCenter();
		//obs->outputIonCenter();
		vars->time+=(logger*dt);
    }
    if (itime%OBSERVE==0) {
        obs->display();
        obs->export_dump_close();
    }
	if ((itime+1)%OBSERVE==0) {
		flags->eflag=true;
		vars->Uzero();
    }
  	verlet();
	}
}
