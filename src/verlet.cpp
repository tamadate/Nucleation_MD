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

	update_position();
	getIonCenterProp();

	thermo->computeZeta(dt);

	check_pairlist();
	vars->totalVirial=0;
	for (auto &a : InterInter) a->compute(vars,flags);
	for (auto &a : IntraInter) a->compute(vars,flags);

	velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))

	thermo->tempControl(dt);
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
		for(int nth=0;nth<Nth;nth++){
			a.fx=a.fy=a.fz=0;
		}
	}
 	 for (auto &i : vars->vapor_in) {
		for (auto &a : vars->vapors[i].inAtoms){
			a.qx += a.px * dt;
			a.qy += a.py * dt;
			a.qz += a.pz * dt;
		  a.fx=a.fy=a.fz=0.0;
			for(int nth=0;nth<Nth;nth++){
				a.fx=a.fy=a.fz=0;
			}
		}
  }

  	for (auto &i : vars->gas_in) {
		for (auto &a : vars->gases[i].inAtoms){
			a.qx += a.px * dt;
			a.qy += a.py * dt;
			a.qz += a.pz * dt;
		  	a.fx=a.fy=a.fz=0.0;
			for(int nth=0;nth<Nth;nth++){
				a.fx=a.fy=a.fz=0;
			}
		}
  }
	vars->times.tpos+=omp_get_wtime();
}


