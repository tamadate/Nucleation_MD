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
MD::run_diff(char** argv) {
/****Thermal relaxation****/
    const int logger=100;
	flags->inter_vi=0;
  flags->inter_vg=0;
	flags->force_lj=0;
	setPotential(flags);
	for (auto &a : IntraInter) {a->printName();}
	cout<<endl;
	for (auto &a : InterInter) {a->printName();}
	for (itime=0; itime < pp->step_relax; itime++) {
		if (itime%OBSERVE==0) {
			display(1);
			export_dump_close();
			vars->time=(itime*dt);
		}
    verlet();
		//if ((itime+1)%OBSERVE==0) flags->eflag=1;
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
	flags->inter_vi=1;
	flags->inter_vg=1;
	flags->force_lj=1;
	flags->inter_vv=1;
	setPotential(flags);
	for (auto &a : vars->ions) a.qx-=ion_r[0], a.qy-=ion_r[1], a.qz-=ion_r[2];
	for (auto &a : vars->gases) a.qx-=ion_r[0], a.qy-=ion_r[1], a.qz-=ion_r[2];
	analysis_ion();
	make_pair();
	margin_length = MARGIN;
	setNVE();

/****Main simulaiton****/
	for (itime=0; itime<Noftimestep; itime++) {
        if(itime%logger==0)   {
            analysis_gas();
            output();
            //output_gas();
            vars->time+=(logger*dt);
            if (itime%OBSERVE==0) {
                display(0);
                export_dump_close();
            }
        }
			if ((itime+1)%OBSERVE==0) flags->eflag=1;
    	verlet();
	}
}
