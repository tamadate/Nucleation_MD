#include "md.hpp"
/*########################################################################################

-----Output-----

#######################################################################################*/

/**********************************initialization******************************************/
void
MD::output_initial(void){
	sprintf(filepath, "ion_%d_%d.dat", int(T), int(calculation_number)); 
	FILE*f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "%d_%d_relax.dat", int(T), int(calculation_number)); 
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "gas_%d_%d.dat", int(T), int(calculation_number)); 
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "vapor_collision_%d.dat", int(calculation_number)); 
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "vapor_in_%d.dat", int(calculation_number)); 
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "vapor_out_%d.dat", int(calculation_number)); 
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "gas_collision_%d.dat", int(calculation_number)); 
	f=fopen(filepath, "w");
	fclose(f);
}


void
MD::output(void){
	sprintf(filepath, "ion_%d_%d.dat", int(T), int(calculation_number)); 
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, ion_r[0], ion_r[1], ion_r[2], ion_v[0], ion_v[1], ion_v[2]);
	fclose(f);
}

void
MD::output_gas(void){
	sprintf(filepath, "gas_%d_%d.dat", int(T), int(calculation_number)); 
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, gas_r[0], gas_r[1], gas_r[2], gas_v[0], gas_v[1], gas_v[2]);
	fclose(f);
}

void
MD::output_gas_collision(long int initime){
	sprintf(filepath, "gas_collision_%d.dat", int(calculation_number)); 
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%e %e\n", initime*dt, itime*dt);
	fclose(f);
}

void
MD::output_vapor_collision(long int initime){
	sprintf(filepath, "vapor_collision_%d.dat", int(calculation_number)); 
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%e %e\n", initime*dt, itime*dt);
	fclose(f);
}

void
MD::Ovin(int i){
	sprintf(filepath, "vapor_in_%d.dat", int(calculation_number)); 
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,itime*dt,vars->vapors[i].qx-ion_r[0],vars->vapors[i].qy-ion_r[1],vars->vapors[i].qz-ion_r[2],vars->vapors[i].px-ion_v[0],vars->vapors[i].py-ion_v[1],vars->vapors[i].pz-ion_v[2]);
	fclose(f);
}

void
MD::Ovout(int i){
	sprintf(filepath, "vapor_out_%d.dat", int(calculation_number)); 
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,itime*dt,vars->vapors[i].qx-ion_r[0],vars->vapors[i].qy-ion_r[1],vars->vapors[i].qz-ion_r[2],vars->vapors[i].px-ion_v[0],vars->vapors[i].py-ion_v[1],vars->vapors[i].pz-ion_v[2]);
	fclose(f);
}

void
MD::output_temp(double gastemp, double iontemp){
	sprintf(filepath, "%d_%d_relax.dat", int(T), int(calculation_number)); 
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%f\t%f\t%f\n", vars->time/1e6, gastemp, iontemp);
	fclose(f);
}


void
MD::display(int output_ONOFF){
    double gastemp = obs->gas_temperature(vars);
    double gastemp_total = obs->gas_total_temperature(vars); 
    double iontemp = obs->ion_temperature(vars);
    double virial=0;//ters->compute_tersoff_virial(vars)/3.0/V*Cpress;
    double gaspress=(kb*Nof_around_gas*gastemp + vars->totalVirial/3.0*6.95e-21)/(V*1e-30);
    double U = vars->totalPotential;
    double K = obs->ion_kinetic_energy(vars);
    std::cout << "----------------------TIME = " << vars->time/1000.0 << " ps-------------------------" << endl;
    printf("Kinetic energy = %1.2e		Temp(gas) = %f\nPotential energy = %1.2e		Temp(ion) = %f\nTotal energy = %1.2e			Temp(gas total) = %f\nPress(gas) = %f\n", K, gastemp, U, iontemp, K+U, gastemp_total, gaspress/101300.0); 
/*    FILE*f=fopen("T-E.dat", "a");
    fprintf(f, "%e\t%e\n", iontemp, U);
    fclose(f);*/
}

void
MD::export_dump(void) {
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	int Nion=int(vars->ions.size());
	int Nvapor=0;
	Nvapor+=int(vars->vapor_out.size());
	Nvapor+=vars->vapors[0].inAtoms.size()*int(vars->vapor_in.size());
	int Ngas=0;
	Ngas+=int(vars->gas_out.size());
	Ngas+=vars->gases[0].inAtoms.size()*int(vars->gas_in.size());

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz\n", count, Ngas+Nion+Nvapor, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {X=ion_r[0];Y=ion_r[1];Z=ion_r[2];}

	int ID=0;
	for (auto &a : vars->ions) {
		fprintf(f,"%d %d %f %f %f %f %f %f\n",ID,a.type,a.qx-X,a.qy-Y,a.qz-Z,a.px,a.py,a.pz);
		ID++;
	}
    for (auto &a : vars->gas_out) {
		fprintf(f, "%d %d %f %f %f %f %f %f\n", ID, 222, vars->gases[a].qx-X, vars->gases[a].qy-Y, vars->gases[a].qz-Z, vars->gases[a].px, vars->gases[a].py, vars->gases[a].pz);
		ID++;
	}
    for (auto &a : vars->gas_in) {
		for (auto &b : vars->gases[a].inAtoms) {
			fprintf(f, "%d %d %f %f %f %f %f %f\n", ID, b.type, b.qx-X, b.qy-Y, b.qz-Z, b.px, b.py, b.pz);
			ID++;
		}
	}
    for (auto &a : vars->vapor_out) {
		fprintf(f, "%d %d %f %f %f %f %f %f\n", ID, 333, vars->vapors[a].qx-X, vars->vapors[a].qy-Y, vars->vapors[a].qz-Z, vars->vapors[a].px, vars->vapors[a].py, vars->vapors[a].pz);
		ID++;
	}
    for (auto &a : vars->vapor_in) {
		for (auto &b : vars->vapors[a].inAtoms) {
			fprintf(f, "%d %d %f %f %f %f %f %f\n", ID, b.type, b.qx-X, b.qy-Y, b.qz-Z, b.px, b.py, b.pz);
			ID++;
		}
	}
	fclose(f);
}


void
MD::export_dump_close(void) {
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	int Nion=int(vars->ions.size());
	int Nvapor=vars->vapors[0].inAtoms.size()*int(vars->vapor_in.size());;
	int Ngas=vars->gases[0].inAtoms.size()*int(vars->gas_in.size());

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz\n", count, Ngas+Nion+Nvapor, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {X=ion_r[0];Y=ion_r[1];Z=ion_r[2];}
	int ID=0;
	for (auto &a : vars->ions) {
		fprintf(f,"%d %d %f %f %f %f %f %f\n",ID,a.type,a.qx-X,a.qy-Y,a.qz-Z,a.px,a.py,a.pz);
		ID++;
	}
    for (auto &a : vars->gas_in) {
		for (auto &b : vars->gases[a].inAtoms) {
			fprintf(f, "%d %d %f %f %f %f %f %f\n", ID, b.type, b.qx-X, b.qy-Y, b.qz-Z, b.px, b.py, b.pz);
			ID++;
		}
	}
    for (auto &a : vars->vapor_in) {
		for (auto &b : vars->vapors[a].inAtoms) {
			fprintf(f, "%d %d %f %f %f %f %f %f\n", ID, b.type, b.qx-X, b.qy-Y, b.qz-Z, b.px, b.py, b.pz);
			ID++;
		}
	}
	fclose(f);
}



