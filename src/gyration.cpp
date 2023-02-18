#include "md.hpp"

void
MD::gyration_out(MD *md2){
	FILE*f=fopen(pp->gyration_path, "a");
	double numerator=0;
	for(auto &a : vars->ions){
		double dx=a.qx-vars->IonX[0];
		double dy=a.qy-vars->IonX[1];
		double dz=a.qz-vars->IonX[2];
		double rsq=dx*dx+dy*dy+dz*dz;
		numerator+=a.mass*rsq;
	}
	fprintf(f, "%f\t%f\t",vars->time,sqrt(numerator/pp->Mion));
	numerator=0;
	for(auto &a : md2->vars->ions){
		double dx=a.qx-md2->vars->IonX[0];
		double dy=a.qy-md2->vars->IonX[1];
		double dz=a.qz-md2->vars->IonX[2];
		double rsq=dx*dx+dy*dy+dz*dz;
		numerator+=a.mass*rsq;
	}
	fprintf(f, "%f\n",sqrt(numerator/md2->pp->Mion));

	fclose(f);
}
