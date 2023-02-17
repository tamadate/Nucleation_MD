#include "potential.hpp"


void
PotentialGasIntra::compute(Variables *vars, FLAG *flags){
	vars->times.tgas-=omp_get_wtime();
	for(auto &i : vars->gas_in){
		Atom *g1= & vars->gases[i].inAtoms[0];
		Atom *g2= & vars->gases[i].inAtoms[1];
    double dx = g1->qx - g2->qx;
    double dy = g1->qy - g2->qy;
    double dz = g1->qz - g2->qz;
    double rsq = (dx * dx + dy * dy + dz * dz);
    double r = sqrt(rsq);
		double dr= r - 1.098;
    double rk = 1221.7 * dr;
    double force_bond_harmonic = -2.0*rk/r;
    g1->fx += force_bond_harmonic * dx;
    g1->fy += force_bond_harmonic * dy;
    g1->fz += force_bond_harmonic * dz;
    g2->fx -= force_bond_harmonic * dx;
    g2->fy -= force_bond_harmonic * dy;
    g2->fz -= force_bond_harmonic * dz;
		if(flags->eflag) vars->Utotal.Ugas+=rk*dr;
	}
	vars->times.tgas+=omp_get_wtime();
}
