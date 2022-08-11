//------------------------------------------------------------------------
#include "potential.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on ion-ion (Coulombic+LJ)
*/
/////////////////////////////////////////////////////////////////////
void
Potential::compute(Variables *vars, FLAG *flags){
	Ion *ions = vars->ions.data();

	for (auto &a : vars->ion_pairs){
		int i=a.i;
		int j=a.j;
		double dx = ions[i].qx - ions[j].qx;
		double dy = ions[i].qy - ions[j].qy;
		double dz = ions[i].qz - ions[j].qz;
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r2inv = 1/rsq;
		int type1=ions[i].type;
		int type2=ions[j].type;
		double r6inv = r2inv * r2inv * r2inv;
		double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);	
		double force_coul = qqrd2e * ions[i].charge * ions[j].charge * sqrt(r2inv);
		double force_pair = (force_lj + force_coul)*r2inv;
		ions[i].fx += force_pair * dx;
		ions[i].fy += force_pair * dy;
		ions[i].fz += force_pair * dz;
		ions[j].fx -= force_pair * dx;
		ions[j].fy -= force_pair * dy;
		ions[j].fz -= force_pair * dz;
		if(flags->eflag) {
			vars->totalPotential+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			vars->totalPotential+=qqrd2e * ions[i].charge * ions[j].charge * sqrt(r2inv);
		}
/*		if(flags->vflag){
			vars->totalVirial+=force_lj+force_coul;
		}*/
	}
}


void
PotentialGasIntra::compute(Variables *vars, FLAG *flags){
	Gas *gases = vars->gases.data();
	for(auto &i : vars->gas_in){
		if(i<Nof_around_gas){
		    int j=i+Nof_around_gas;
		    double dx = gases[i].qx - gases[j].qx;
		    double dy = gases[i].qy - gases[j].qy;
		    double dz = gases[i].qz - gases[j].qz;
		    double rsq = (dx * dx + dy * dy + dz * dz);
		    double r = sqrt(rsq);
			double dr= r - 1.098;
		    double rk = 1221.7 * dr;
		    double force_bond_harmonic = -2.0*rk/r;
		    gases[i].fx += force_bond_harmonic * dx;
		    gases[i].fy += force_bond_harmonic * dy;
		    gases[i].fz += force_bond_harmonic * dz;
		    gases[j].fx -= force_bond_harmonic * dx;
		    gases[j].fy -= force_bond_harmonic * dy;
		    gases[j].fz -= force_bond_harmonic * dz;
			if(flags->eflag) vars->totalPotential+=rk*dr;

		}
	}
}




