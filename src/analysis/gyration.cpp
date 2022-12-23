#include "../md.hpp"

void
MD::gyration_out(MD *md2){
	FILE*f=fopen(pp->gyration_path, "a");
	double numerator=0;
	for(auto &a : vars->Molecules[0].inAtoms){
		double dx=a.qx-vars->Molecules[0].qx;
		double dy=a.qy-vars->Molecules[0].qy;
		double dz=a.qz-vars->Molecules[0].qz;
		double rsq=dx*dx+dy*dy+dz*dz;
		numerator+=a.mass*rsq;
	}
	fprintf(f, "%f\t%f\t",vars->time,sqrt(numerator/pp->Mion));
	fclose(f);
}
