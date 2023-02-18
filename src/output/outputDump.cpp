#include "../md.hpp"

void
Observer::outputDump(void) {
	int count = vars->time;
	FILE*f=fopen(fileDump, "a");
	int Nion=int(vars->ions.size());
	int Nvapor=0;
	Nvapor+=int(vars->vapor_out.size());
	Nvapor+=vars->vapors[0].inAtoms.size()*int(vars->vapor_in.size());
	int Ngas=0;
	Ngas+=int(vars->gas_out.size());
	Ngas+=vars->gases[0].inAtoms.size()*int(vars->gas_in.size());

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz\n", count, Ngas+Nion+Nvapor, 
        -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(dump_fix==1) {X=vars->ion_r[0];Y=vars->ion_r[1];Z=vars->ion_r[2];}

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
Observer::outputDumpClose(void) {
	int count = vars->time;
	FILE*f=fopen(fileDump, "a");
	int Nion=vars->ions.size();
	int Nvapor=vars->vapors[0].inAtoms.size()*int(vars->vapor_in.size());;
	int Ngas=vars->gases[0].inAtoms.size()*int(vars->gas_in.size());

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz\n", 
        count, Ngas+Nion+Nvapor, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(dump_fix==1) {X=vars->ion_r[0];Y=vars->ion_r[1];Z=vars->ion_r[2];}
	int ID=0;
	for (auto &a : vars->ions) {
		fprintf(f,"%d %s %f %f %f %f %f %f\n",ID,vars->atypes[a.type].name.c_str(),a.qx-X,a.qy-Y,a.qz-Z,a.px,a.py,a.pz);
		ID++;
	}
    for (auto &a : vars->gas_in) {
		for (auto &b : vars->gases[a].inAtoms) {
			fprintf(f, "%d %s %f %f %f %f %f %f\n", ID, (vars->atypes[0].name+"(g)").c_str(), b.qx-X, b.qy-Y, b.qz-Z, b.px, b.py, b.pz);
			ID++;
		}
	}
    for (auto &a : vars->vapor_in) {
		for (auto &b : vars->vapors[a].inAtoms) {
			fprintf(f, "%d %s %f %f %f %f %f %f\n", ID, (vars->atypes[b.type].name+"(v)").c_str(), b.qx-X, b.qy-Y, b.qz-Z, b.px, b.py, b.pz);
			ID++;
		}
	}
	fclose(f);
}
