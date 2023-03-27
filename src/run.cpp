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
	// Thermal relaxation
	for (itime=0; itime < con->step_relax; itime++) {
		vars->time=itime*con->dt;
		if (itime%obs->OBSERVE==0) {
			vars->eflag=true;
			vars->Uzero();
		}
		
		verlet();

		if (itime%obs->OBSERVE==0) {
			obs->display();
			obs->outputDumpClose();
			vars->eflag=false;
		}
	}


	// Delete thermostat and set NVE 
	int ifunc=0;
	for(auto &a : funcs){
		if(a -> type == 1){
			funcs.erase(funcs.begin()+ifunc);
		}
		else{
			ifunc++;
		}
	}
	funcs.push_back(new NVE());

	// Main simulaiton
	for (itime=0; itime<con->Noftimestep; itime++) {
		vars->time=itime*con->dt;
		if (itime%obs->OBSERVE==0) {
			vars->eflag=true;
			vars->Uzero();
		}

		verlet();

		if (itime%obs->OBSERVE==0) {
			obs->display();
			obs->outputDumpClose();
			vars->eflag=false;
		}
	}
}
