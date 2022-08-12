//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*	
	- Compute a domain (NVE or NVT)
	- v(t) -> v(t+dt/2)
	- Calculate velocity and potition of ion's center of mass.
	- Temperature control (velocity scaling, this case is NVT)
	- r(t) -> r(t+dt)
	- apply periodic boundary condition
	- Determine whethere update pair list or not
	- Calculate ion intraatmic interaction
	- Calculate ion-gas interatominc interaction 
	- v(t+dt/2) -> v(t)
*/
/////////////////////////////////////////////////////////////////////
void
MD::verlet(void) {
	velocity_calculation(); //	v(t) -> v(t+dt/2) using F(x(t))
	analysis_ion();
	if(flags->fix_cell_center==1) fix_cell_center();
	if(flags->velocity_scaling==1)	velocity_scaling();
    
	update_position();	

	if(flags->nose_hoover_ion==1)	nosehoover_zeta();
	if(flags->nose_hoover_gas==1)	nosehoover_zeta_gas();
	check_pairlist();
	vars->totalPotential=0;
	vars->totalVirial=0;

	for (auto &a : InterInter) a->compute(vars,flags);
	for (auto &a : IntraInter) a->compute(vars,flags);

	velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))
	if(flags->nose_hoover_ion==1)	nosehoover_ion();
	if(flags->nose_hoover_gas==1)	nosehoover_gas();
}

/////////////////////////////////////////////////////////////////////
/*	
	- Update velocity (half of a time step, dt/2)
*/
/////////////////////////////////////////////////////////////////////
void 
MD::velocity_calculation(void) {
	Gas *gases = vars->gases.data();
	double const Coeff=0.5*dt*4.184e-4;
	for (auto &a : vars->ions) {
		double Coeff2=Coeff/a.mass;
		a.px += a.fx *Coeff2;
	    a.py += a.fy *Coeff2;
	    a.pz += a.fz *Coeff2;
	}
    for (auto &i : vars->vapor_in) {
		for (auto &a : vars->vapors[i].inAtoms){
			double Coeff2=Coeff/a.mass;
			a.px += a.fx * Coeff2;
			a.py += a.fy * Coeff2;
			a.pz += a.fz * Coeff2;
		}
    }
    for (auto &i : vars->gas_in){
		double Coeff2=Coeff/gases[i].mass;
        gases[i].px += gases[i].fx *Coeff2;
        gases[i].py += gases[i].fy *Coeff2;
        gases[i].pz += gases[i].fz *Coeff2;
	}
}

/////////////////////////////////////////////////////////////////////
/*	
	- Update velocity (a time step, dt)
*/
/////////////////////////////////////////////////////////////////////
void
MD::update_position(void) {
	for (auto &a : vars->ions) {
		a.qx += a.px * dt;
		a.qy += a.py * dt;
		a.qz += a.pz * dt;
		a.fx=a.fy=a.fz=0.0;
	}
    for (auto &i : vars->vapor_in) {
		for (auto &a : vars->vapors[i].inAtoms){
			a.qx += a.px * dt;
			a.qy += a.py * dt;
			a.qz += a.pz * dt;
		    a.fx=a.fy=a.fz=0.0;    
		}
    }

	Gas *gases = vars->gases.data();
    for (auto &i : vars->gas_in) {
        gases[i].qx += gases[i].px * dt;
        gases[i].qy += gases[i].py * dt;
        gases[i].qz += gases[i].pz * dt;
        gases[i].fx=gases[i].fy=gases[i].fz=0.0;    
    }
}


/////////////////////////////////////////////////////////////////////
/*	
	- Fix center of domain (v-v_center)
*/
/////////////////////////////////////////////////////////////////////
void
MD::fix_cell_center(void) {
	for (auto &a : vars->ions) {
		a.px -= ion_v[0];
		a.py -= ion_v[1];
		a.pz -= ion_v[2];
	}
}
