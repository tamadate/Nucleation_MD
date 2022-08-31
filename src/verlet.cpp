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
	//if(flags->nose_hoover_gas==1)	nosehoover_zeta_gas();
	check_pairlist();
	vars->totalVirial=0;
	for (auto &a : InterInter) a->compute(vars,flags);
	for (auto &a : IntraInter) a->compute(vars,flags);

	velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))

	if(flags->nose_hoover_ion==1)	nosehoover_ion();
	//(flags->nose_hoover_gas==1)	nosehoover_gas();
	flags->eflag=0;
}

/////////////////////////////////////////////////////////////////////
/*
	- Update velocity (half of a time step, dt/2)
*/
/////////////////////////////////////////////////////////////////////
void
MD::velocity_calculation(void) {
	vars->times.tvel-=omp_get_wtime();
	Molecule *gases = vars->gases.data();
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
		for (auto &a : vars->gases[i].inAtoms){
			double Coeff2=Coeff/a.mass;
			a.px += a.fx * Coeff2;
			a.py += a.fy * Coeff2;
			a.pz += a.fz * Coeff2;
		}
	}
	vars->times.tvel+=omp_get_wtime();
}

/////////////////////////////////////////////////////////////////////
/*
	- Update velocity (a time step, dt)
*/
/////////////////////////////////////////////////////////////////////
void
MD::update_position(void) {
	vars->times.tpos-=omp_get_wtime();
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

  for (auto &i : vars->gas_in) {
		for (auto &a : vars->gases[i].inAtoms){
			a.qx += a.px * dt;
			a.qy += a.py * dt;
			a.qz += a.pz * dt;
		  a.fx=a.fy=a.fz=0.0;
		}
  }
	vars->times.tpos+=omp_get_wtime();
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
