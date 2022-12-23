#include "../md.hpp"
/*########################################################################################

-----Output-----

#######################################################################################*/

/**********************************initialization******************************************/
void
MD::output_initial(void){
	sprintf(filepath, "ion_%d_%d.dat", int(T), int(calculation_number));
	FILE*f=fopen(filepath, "w");
	fclose(f);
	/*sprintf(filepath, "%d_%d_relax.dat", int(T), int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "gas_%d_%d.dat", int(T), int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);*/
	sprintf(filepath, "vapor_collision_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "vapor_in_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "vapor_out_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	/*sprintf(filepath, "gas_collision_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);*/
	sprintf(filepath, "K_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "U_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
}


void
MD::output(void){
	sprintf(filepath, "ion_%d_%d.dat", int(T), int(calculation_number));
	FILE*f=fopen(filepath, "a");
	Molecule *ions = vars->Molecules.data();
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, ions[0].qx, ions[0].qy, ions[0].qz, ions[0].px, ions[0].py, ions[0].pz);
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
	Molecule *mols=vars->Molecules.data();
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,itime*dt,mols[i].qx-mols[0].qx,mols[i].qy-mols[0].qy,mols[i].qz-mols[0].qz,\
	mols[i].px-mols[0].px,mols[i].py-mols[0].py,mols[i].pz-mols[0].pz);
	fclose(f);
}

void
MD::Ovout(int i){
	sprintf(filepath, "vapor_out_%d.dat", int(calculation_number));
	FILE*f=fopen(filepath, "a");
	Molecule *mols=vars->Molecules.data();
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,itime*dt,mols[i].qx-mols[0].qx,mols[i].qy-mols[0].qy,mols[i].qz-mols[0].qz,\
	mols[i].px-mols[0].px,mols[i].py-mols[0].py,mols[i].pz-mols[0].pz);
	fclose(f);
}

void
MD::display(int output_ONOFF){
		obs->computeProps(vars,0);
		obs->computeProps(vars,1);
		obs->computeProps(vars,2);
    double virial=0;//ters->compute_tersoff_virial(vars)/3.0/V*Cpress;
    double gaspress=(kb*Nof_around_gas*obs->Tout[1] + vars->totalVirial/3.0*6.95e-21)/(V*1e-30);
    double U = vars->Usum();
		double Kin=0;
		for(auto &kin : obs->Kin) Kin+=kin;
		double Kout=0;
		for(auto &kout : obs->Kout) Kout+=kout;

    std::cout << "-------------------TIME = " << vars->time/1000.0 << " ps----------------------" << endl;
		cout<<"*Inside propeties"<<endl;
		printf("  %-15s%-12s%-12s%-12s\n", "Prop", "Ion", "Gas", "Vapor");
		printf("  %-15s%-12.2e%-12.2e%-12.2e\n", "K-Energy", obs->Kin[0], obs->Kin[1], obs->Kin[2]);
		printf("  %-15s%-12.2f%-12.2f%-12.2f\n", "Temperature", obs->Tin[0], obs->Tin[1], obs->Tin[2]);
		printf("  %-15s%-12.2e%-12.2e%-12.2e\n", "Uintra", vars->Utotal.Uion, vars->Utotal.Ugas, vars->Utotal.Uvap);
		printf("  %-15s%-12.2e%-12.2e%-12.2e\n", "Uinter(Ion)", 0.0, vars->Utotal.Ugi, vars->Utotal.Uvi);
		printf("  %-15s%-12s%-12.2e%-12.2e\n", "Uinter(Gas)", "-", vars->Utotal.Ugg, vars->Utotal.Uvg);
		printf("  %-15s%-12s%-12s%-12.2e\n", "Uinter(Vapor)", "-", "-", vars->Utotal.Uvv);

		cout<<"*Out side propeties"<<endl;
		printf("  %-15s%-12s%-12s%-12s\n", "Prop", "Ion", "Gas", "Vapor");
		printf("  %-15s%-12s%-12.2e%-12.2e\n", "K-Energy", "-", obs->Kout[1], obs->Kout[2]);
		printf("  %-15s%-12s%-12.2f%-12.2f\n", "Temperature", "-", obs->Tout[1], obs->Tout[2]);

		cout<<"*System total propeties"<<endl;
		printf("  K = %1.2e	U = %1.2e	Press = %f\n",Kin+Kout, U, gaspress/101300.0);

		cout<<"*Times"<<endl;
		printf("  tion  = %1.1f s	tgas = %1.1f s 	tvap = %1.1f s\n",vars->times.tion,vars->times.tgas,vars->times.tvap);
		printf("  tvi   = %1.1f s	tgi  = %1.1f s	tvg  = %1.1f s	tvv  = %1.1f s\n",vars->times.tvi,vars->times.tgi,vars->times.tvg,vars->times.tvv);
		printf("  tpair = %1.1f s	tpos = %1.1f s	tvel = %1.1f s	tetc = %1.1f s\n",vars->times.tpair,vars->times.tpos,vars->times.tvel,vars->times.tetc);
		printf("  tpot  = %1.1f s	ttot = %1.1f s\n",(vars->times.tvi+vars->times.tgi+vars->times.tvv+vars->times.tvg+vars->times.tion+vars->times.tgas+vars->times.tvap), omp_get_wtime()-startTime);
		printf("  NCPU = %d\n",Nth);

		cout <<endl;

		sprintf(filepath, "K_%d.dat", int(calculation_number));
		FILE*f=fopen(filepath, "a");
		fprintf(f,"%e %e %e %e %e %e\n",vars->time,obs->Kin[0],obs->Kin[1],obs->Kin[2],obs->Kout[1],obs->Kout[2]);
		fclose(f);

		sprintf(filepath, "U_%d.dat", int(calculation_number));
		f=fopen(filepath, "a");
		fprintf(f,"%e %e %e %e %e %e %e %e %e\n",vars->time,vars->Utotal.Uion,vars->Utotal.Ugas,vars->Utotal.Uvap,vars->Utotal.Ugi,vars->Utotal.Ugg,vars->Utotal.Uvg,vars->Utotal.Uvi,vars->Utotal.Uvv);
		fclose(f);
}

void
MD::exportDumpIn(void) {
	Molecule *mols = vars->Molecules.data();
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	int Natom=0;
	for (auto &a : vars->Molecules){
		Natom+=a.inAtoms.size();
	}

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz\n",
	count, Natom, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {
		X=vars->Molecules[0].qx;
		Y=vars->Molecules[0].qy;
		Z=vars->Molecules[0].qz;
	}

	int ID=0;
	for (auto &a : vars->Molecules[0].inAtoms) {
		fprintf(f,"%d %s %f %f %f %f %f %f\n",ID,vars->atypes[a.type].name.c_str(),a.qx-X,a.qy-Y,a.qz-Z,a.px,a.py,a.pz);
		ID++;
	}
	for (auto i : vars->MolID[1]) {
		for (auto &b : mols[i].inAtoms) {
			fprintf(f, "%d %s %f %f %f %f %f %f\n", ID, (vars->atypes[b.type].name+"(g)").c_str(), b.qx-X, b.qy-Y, b.qz-Z, b.px, b.py, b.pz);
			ID++;
		}
	}
	for (auto i : vars->MolID[2]) {
		for (auto &b : mols[i].inAtoms) {
			fprintf(f, "%d %s %f %f %f %f %f %f\n", ID, (vars->atypes[b.type].name+"(v)").c_str(), b.qx-X, b.qy-Y, b.qz-Z, b.px, b.py, b.pz);
			ID++;
		}
	}
	fclose(f);
}

// secret command
void
MD::exportDumpOut(void) {
	Molecule *mols = vars->Molecules.data();
	int count = vars->time;
	FILE*f=fopen("out.dump", "a");
	int Natom=vars->Molecules.size();

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz\n",
	count, Natom, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {
		X=vars->Molecules[0].qx;
		Y=vars->Molecules[0].qy;
		Z=vars->Molecules[0].qz;
	}

	int ID=0;
	for(int i=3;i<3;i++){
		for (auto j : vars->MolID[i]) {
			fprintf(f,"%d %d %f %f %f %f %f %f\n",ID,i,mols[j].qx-X,mols[j].qy-Y,mols[j].qz-Z,mols[j].px,mols[j].py,mols[j].pz);
			ID++;
		}
	}
	fclose(f);
}
