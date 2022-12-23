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
	if(flags->velocity_scaling==1)	velocity_scaling();

	update_position();

	if(flags->nose_hoover_ion==1)	nosehoover_zeta();
	check_pairlist();
	vars->totalVirial=0;
	for (auto &a : InterInter) a->compute(vars,flags);
	for (auto &a : IntraInter) a->compute(vars,flags);

	//forceCombine();

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
	double const Coeff=0.5*dt*4.184e-4;
	int i=0;
	for (auto &mol : vars->Molecules){
		if(vars->Region[i]==CG){
			double Coeff2=Coeff/mol.mass;
			mol.px += mol.fx *Coeff2;
			mol.py += mol.fy *Coeff2;
			mol.pz += mol.fz *Coeff2;
		}
		else{
			for (auto &a : mol.inAtoms){
				double Coeff2=Coeff/a.mass;
				a.px += a.fx *Coeff2;
				a.py += a.fy *Coeff2;
				a.pz += a.fz *Coeff2;
			}
		}
		i++;
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
	int i=0;
	for (auto &mol : vars->Molecules){
		if(vars->Region[i]==CG){
			mol.qx += mol.px * dt;
			mol.qy += mol.py * dt;
			mol.qz += mol.pz * dt;
			mol.fx=mol.fy=mol.fz=0.0;
		}
		else{
			for (auto &a : mol.inAtoms){
				a.qx += a.px * dt;
				a.qy += a.py * dt;
				a.qz += a.pz * dt;
				a.fx=a.fy=a.fz=0.0;
			}
		}
		i++;
	}
	vars->times.tpos+=omp_get_wtime();

	//boundary_scaling_ion_move();
	updateInCenters();
	boundary_scaling_gas_move();
	boundary_scaling_vapor_move();
}
