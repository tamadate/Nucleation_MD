#include "md.hpp"

void
MD::run(char** argv) {
	// Thermal relaxation
	for (itime=0; itime < con->step_relax; itime++) {
		vars->time=itime*con->dt;	// time = (iteration) * (time step, dt)
		// observer?
		if (itime%obs->OBSERVE==0) {
			vars->eflag=true;
			vars->Uzero();
		}
		
		verlet(); // one step forward iteration

		// observer?
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

		verlet();	// one step forward iteration

		if (itime%obs->OBSERVE==0) {
			obs->display();
			obs->outputDumpClose();
			vars->eflag=false;
		}
	}
}
