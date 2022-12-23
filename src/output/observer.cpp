//------------------------------------------------------------------------
#include "observer.hpp"
//------------------------------------------------------------------------
void
Observer::computeProps(Variables *vars,int molID){
	Kin[molID]=0;
	double Nin=0;
	Molecule *mols = vars -> Molecules.data();

	for(auto i : vars->MolID[molID]){
		for (auto &at : mols[i].inAtoms){
			Kin[molID] += at.px * at.px * at.mass;
			Kin[molID] += at.py * at.py * at.mass;
			Kin[molID] += at.pz * at.pz * at.mass;
			Nin++;
		}
	}
	Kin[molID]*= (0.5 * real_to_kcalmol);
	Tin[molID]=Kin[molID]/double(Nin)*coeff;

	Kout[molID]=0;
	double Nout=0;

	for(auto i : vars->MolID[molID]){
		Kout[molID] += mols[i].px * mols[i].px * mols[i].mass;
		Kout[molID] += mols[i].py * mols[i].py * mols[i].mass;
		Kout[molID] += mols[i].pz * mols[i].pz * mols[i].mass;
		Nout++;
	}
	Kout[molID]*= (0.5 * real_to_kcalmol);
	Tout[molID]=Kout[molID]/double(Nout)*coeff;

};



double
Observer::pressure(Variables *vars, std::vector<Pair> &pairs, double Treal, double virial,double p,double T) {
	double phi = 0.0;
	/*const int ps = pairs.size();
	Gas *gases = vars->gases.data();
	for (int k = 0; k < ps; k++) {
		const int i = pairs[k].i;
		const int j = pairs[k].j;
		double dx = gases[j].qx - gases[i].qx;
		double dy = gases[j].qy - gases[i].qy;
		double dz = gases[j].qz - gases[i].qz;
		adjust_periodic(dx, dy, dz);
		double r2 = (dx * dx + dy * dy + dz * dz);
		double r2inv= 1/r2;
		int type1=gases[i].type;
		int type2=gases[j].type;
		if (r2 < CL2){
			double r6inv = r2inv * r2inv * r2inv;
			phi += r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
		}
	}*/
	phi = phi * Cpress;
	return  p/T*Treal + phi;
}
