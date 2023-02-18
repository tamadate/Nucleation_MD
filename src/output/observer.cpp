//------------------------------------------------------------------------
#include "observer.hpp"
//------------------------------------------------------------------------


Observer::Observer(Variables *VARS, int calculation_number){
	vars=VARS;
	dump_fix=true;
	startTime=omp_get_wtime();
	sprintf(fileIonCenter, "ion_%d.dat", int(calculation_number));
	FILE*f=fopen(fileIonCenter, "w");
	fclose(f);
	sprintf(fileVaporIn, "vapor_in_%d.dat", int(calculation_number));
	f=fopen(fileVaporIn, "w");
	fclose(f);
	sprintf(fileVaporOut, "vapor_out_%d.dat", int(calculation_number));
	f=fopen(fileVaporOut, "w");
	fclose(f);
	sprintf(fileKinetic, "K_%d.dat", int(calculation_number));
	f=fopen(fileKinetic, "w");
	fclose(f);
	sprintf(filePotential, "U_%d.dat", int(calculation_number));
	f=fopen(filePotential, "w");
	fclose(f);
}


double
Observer::pressure(std::vector<Pair> &pairs, double Treal, double virial,double p,double T) {
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
