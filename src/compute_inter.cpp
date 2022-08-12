//------------------------------------------------------------------------
#include "potential.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on ion-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialGasIon::compute(Variables *vars, FLAG *flags) {
	Atom *gases = vars->gases.data();
	Atom *ions = vars->ions.data();
	const int gs = vars->gas_in.size();
	const int is = vars->ions.size();
	for(auto &a : vars->pairs_gi){
		int i=a.i;
		int j=a.j;
		double dx = gases[i].qx - ions[j].qx;
		double dy = gases[i].qy - ions[j].qy;
		double dz = gases[i].qz - ions[j].qz;
		adjust_periodic(dx, dy, dz);
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r2inv = 1/rsq;
		int type1=gases[i].type;
		int type2=ions[j].type;
		if (rsq < CL2) {
			double r6inv = r2inv * r2inv * r2inv;
			double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
			double force_pair = (force_lj)*r2inv;
			gases[i].fx += force_pair * dx;
			gases[i].fy += force_pair * dy;
			gases[i].fz += force_pair * dz;
			ions[j].fx -= force_pair * dx;
			ions[j].fy -= force_pair * dy;
			ions[j].fz -= force_pair * dz;
			if(flags->eflag) vars->totalPotential+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			vars->totalVirial+=force_lj;
		}	
	}
}


/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on vapor-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialVaporGas::compute(Variables *vars, FLAG *flags) {
	Atom *gases = vars->gases.data();
	Molecule *vapors = vars->vapors.data();

	for(auto I : vars->vapor_in){
		for (auto &av1 : vapors[I].inAtoms){
			for(auto j : vars->gas_in){
				double dx = av1.qx - gases[j].qx;
				double dy = av1.qy - gases[j].qy;
				double dz = av1.qz - gases[j].qz;
				double rsq = (dx * dx + dy * dy + dz * dz);
				double r2inv = 1/rsq;
				int type1=av1.type;
				int type2=gases[j].type;
				double r6inv = r2inv * r2inv * r2inv;
				double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
				double force_pair = (force_lj)*r2inv;
				av1.fx += force_pair * dx;
				av1.fy += force_pair * dy;
				av1.fz += force_pair * dz;
				gases[j].fx -= force_pair * dx;
				gases[j].fy -= force_pair * dy;
				gases[j].fz -= force_pair * dz;
				//	if(flags->eflag) vars->totalPotential+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
				//	vars->totalVirial+=force_lj;
			}
		}	
	}

}


/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on ion-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialVaporIon::compute(Variables *vars, FLAG *flags) {
	Molecule *vapors = vars->vapors.data();
	Atom *ions = vars->ions.data();
	const int is = vars->ions.size();
	for(auto &I : vars->vapor_in){
		for (auto &av : vapors[I].inAtoms){
			for(auto &ai : vars->ions){
				double dx = av.qx - ai.qx;
				double dy = av.qy - ai.qy;
				double dz = av.qz - ai.qz;
				adjust_periodic(dx, dy, dz);
				double rsq = (dx * dx + dy * dy + dz * dz);
				double r2inv = 1/rsq;
				int type1=av.type;
				int type2=ai.type;
				double r6inv = r2inv * r2inv * r2inv;
				double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
				double force_coul = qqrd2e * av.charge * ai.charge * sqrt(r2inv);
				double force_pair = (force_lj + force_coul)*r2inv;
				av.fx += force_pair * dx;
				av.fy += force_pair * dy;
				av.fz += force_pair * dz;
				ai.fx -= force_pair * dx;
				ai.fy -= force_pair * dy;
				ai.fz -= force_pair * dz;
				//	if(flags->eflag) vars->totalPotential+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
				//	vars->totalVirial+=force_lj;
			}
		}	
	}
}

/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on ion-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialVaporVapor::compute(Variables *vars, FLAG *flags) {
	Molecule *vapors = vars->vapors.data();
	const int vs = vars->vapor_in.size();
	if(vs>1){
		for(int i1=0; i1<vs-1; i1++){
			int I=vars->vapor_in[i1];
			for (auto &av1 : vapors[I].inAtoms){
				for(int i2=i1+1; i2<vs; i2++){
					int J=vars->vapor_in[i2];
					for (auto &av2 : vapors[J].inAtoms){
						double dx = av1.qx - av2.qx;
						double dy = av1.qy - av2.qy;
						double dz = av1.qz - av2.qz;
						double rsq = (dx * dx + dy * dy + dz * dz);
						double r2inv = 1/rsq;
						int type1=av1.type;
						int type2=av2.type;
						double r6inv = r2inv * r2inv * r2inv;
						double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
						double force_coul = qqrd2e * av1.charge * av2.charge * sqrt(r2inv);
						double force_pair = (force_lj + force_coul)*r2inv;
						av1.fx += force_pair * dx;
						av1.fy += force_pair * dy;
						av1.fz += force_pair * dz;
						av2.fx -= force_pair * dx;
						av2.fy -= force_pair * dy;
						av2.fz -= force_pair * dz;
						//	if(flags->eflag) vars->totalPotential+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
						//	vars->totalVirial+=force_lj;
					}
				}
			}	
		}
	}
}



/////////////////////////////////////////////////////////////////////
/*	
	- Calculate the ion induced dipole force working between ion and gas 
	molecules. Because ion have a charge, ion induce the gas dipole.
*/
/////////////////////////////////////////////////////////////////////


// charge devided by number of atoms in ions
void
PotentialIonDipole::compute(Variables *vars, FLAG *flags) {
	Atom *gases = vars->gases.data();
	const int gs = vars->gas_in.size();
	const int is = vars->ions.size();
    
    for (int k = 0; k < gs; k++) {
        for (auto &a : vars->ions) {
            int i = vars->gas_in[k];
            double dx = gases[i].qx - a.qx;
            double dy = gases[i].qy - a.qy;
            double dz = gases[i].qz - a.qz;
			adjust_periodic(dx, dy, dz);
			double rsq = (dx * dx + dy * dy + dz * dz);
            
            double r2inv = 1/rsq;
            double ion_dipole = -2*alphagas*r2inv*r2inv*qqrd2e/is/is*zion*zion;
            double force_pair=ion_dipole*r2inv;
            gases[i].fx += force_pair * dx;
            gases[i].fy += force_pair * dy;
            gases[i].fz += force_pair * dz;
            a.fx -= force_pair * dx;
            a.fy -= force_pair * dy;
            a.fz -= force_pair * dz;
			if(flags->eflag) vars->totalPotential += -0.5*alphagas*r2inv*r2inv*qqrd2e/is/is;
		}
	}
}


/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on gas-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialGasGas::compute(Variables *vars, FLAG *flags) {
	Atom *gases = vars->gases.data();
	for(auto &a : vars->pairs_gg){
		int i=a.i;
		int j=a.j;
		double dx = gases[i].qx - gases[j].qx;
		double dy = gases[i].qy - gases[j].qy;
		double dz = gases[i].qz - gases[j].qz;
		adjust_periodic(dx, dy, dz);
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r2inv = 1/rsq;
		int type1=gases[i].type;
		int type2=gases[j].type;
		if (rsq < CL2) {
			double r6inv = r2inv * r2inv * r2inv;
			double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
			double force_pair = (force_lj)*r2inv;
			gases[i].fx += force_pair * dx;
			gases[i].fy += force_pair * dy;
			gases[i].fz += force_pair * dz;
			gases[j].fx -= force_pair * dx;
			gases[j].fy -= force_pair * dy;
			gases[j].fz -= force_pair * dz;
		}	
	}
}

/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on gas-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialEfield::compute(Variables *vars, FLAG *flags) {
    for (auto &a : vars->ions) {
		a.fx+=6.2665e-5*a.charge*Ecoeff[0];
		a.fy+=6.2665e-5*a.charge*Ecoeff[1];
		a.fz+=6.2665e-5*a.charge*Ecoeff[2];
	}
	
}



