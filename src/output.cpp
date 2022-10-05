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
	Molecule *ions = vars->effectiveIn[0].data();
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
	Molecule *vap=vars->effectiveIn[2].data();
	Molecule *ion=vars->effectiveIn[0].data();
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,itime*dt,vap[i].qx-ion[0].qx,vap[i].qy-ion[0].qy,vap[i].qz-ion[0].qz,\
	vap[i].px-ion[0].px,vap[i].py-ion[0].py,vap[i].pz-ion[0].pz);
	fclose(f);
}

void
MD::Ovout(int i){
	sprintf(filepath, "vapor_out_%d.dat", int(calculation_number));
	FILE*f=fopen(filepath, "a");
	Molecule *vap=vars->effectiveIn[2].data();
	Molecule *ion=vars->effectiveIn[0].data();
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,itime*dt,vap[i].qx-ion[0].qx,vap[i].qy-ion[0].qy,vap[i].qz-ion[0].qz,\
	vap[i].px-ion[0].px,vap[i].py-ion[0].py,vap[i].pz-ion[0].pz);
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

		vars->Ucombine();
    std::cout << "----------------------TIME = " << vars->time/1000.0 << " ps-------------------------" << endl;
		cout<<"Inside propeties"<<endl;
    printf("  Kion = %1.2e  Kgas = %1.2e  Kvap = %1.2e\n  Tion = %1.2f  Tgas = %1.2f  Tvap = %1.2f\n\
		  Uion = %1.2e  Ugas = %1.2e  Uvap = %1.2e\n  Ugi = %1.2e  Ugg = %1.2e  Uvi = %1.2e	\n  Uvg = %1.2e	Uvv= %1.2e  \n",\
		obs->Kin[0], obs->Kin[1], obs->Kin[2], obs->Tin[0], obs->Tin[1], obs->Tin[2], vars->Utotal.Uion, vars->Utotal.Ugas,
		vars->Utotal.Uvap, vars->Utotal.Ugi, vars->Utotal.Ugg,	vars->Utotal.Uvi, vars->Utotal.Uvg, vars->Utotal.Uvv);
		cout<<"Out side propeties"<<endl;
		printf("  Kgas = %1.2e  Tgas = %1.2f  Ugas = %1.2e	\n  Kvap = %1.2e  Tvap = %1.2f  Uvap = %1.2e	\n  Kout = %1.2e    Uout = %1.2e	\n",
		obs->Kout[1], obs->Tout[1], 0.0, obs->Kout[2], obs->Tout[2], 0.0, Kout, 0.0);
		cout<<"System propeties"<<endl;
		printf("  K = %1.2e	U = %1.2e	Press = %f\n",Kin+Kout, U, gaspress/101300.0);
		cout<<"Times"<<endl;
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
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	int Natom=0;
	for (auto &a : vars->effectiveIn){
		for (auto &b : a){
			Natom+=b.inAtoms.size();
		}
	}

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz\n",
	count, Natom, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {
		X=vars->effectiveIn[0][0].qx;
		Y=vars->effectiveIn[0][0].qy;
		Z=vars->effectiveIn[0][0].qz;
	}

	int ID=0;
	for (auto &a : vars->effectiveIn[0][0].inAtoms) {
		fprintf(f,"%d %s %f %f %f %f %f %f\n",ID,vars->atypes[a.type].name.c_str(),a.qx-X,a.qy-Y,a.qz-Z,a.px,a.py,a.pz);
		ID++;
	}
	for (auto &a : vars->effectiveIn[1]) {
		for (auto &b : a.inAtoms) {
			fprintf(f, "%d %s %f %f %f %f %f %f\n", ID, (vars->atypes_g.name+"(g)").c_str(), b.qx-X, b.qy-Y, b.qz-Z, b.px, b.py, b.pz);
			ID++;
		}
	}
	for (auto &a : vars->effectiveIn[2]) {
		for (auto &b : a.inAtoms) {
			fprintf(f, "%d %s %f %f %f %f %f %f\n", ID, (vars->atypes_v[b.type].name+"(v)").c_str(), b.qx-X, b.qy-Y, b.qz-Z, b.px, b.py, b.pz);
			ID++;
		}
	}
	fclose(f);
}

// secret command
void
MD::exportDumpOut(void) {
	int count = vars->time;
	FILE*f=fopen("out.dump", "a");
	int Natom=0;
	for (auto &a : vars->effectiveOut) Natom+=a.size();

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz\n",
	count, Natom, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {
		X=vars->effectiveIn[0][0].qx;
		Y=vars->effectiveIn[0][0].qy;
		Z=vars->effectiveIn[0][0].qz;
	}

	int ID=0;
	for(int i=0;i<3;i++){
		for (auto &a : vars->effectiveOut[i]) {
			fprintf(f,"%d %d %f %f %f %f %f %f\n",ID,a.inFlag,a.qx-X,a.qy-Y,a.qz-Z,a.px,a.py,a.pz);
			ID++;
		}
	}
	fclose(f);
}
