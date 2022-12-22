//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*
	make pair list
	- always gas-ion pair list was calculated
	- if gas-gas interaction is ON, make pair_gasgas
	- get max velocity of gas molecule vmax
	- set update loop = margine_length/(vmax*dt)
*/
/////////////////////////////////////////////////////////////////////
void
MD::make_pair(void){
	updateInCenters();
	updateInOut(1);
	updateInOut(2);
	make_pair_gasion();
	make_pair_gasvapor();
	make_pair_vaporvapor();
	make_pair_gasgas();
	// set number of steps for next pair list update
	double vmax2 = 0.0;
	for (auto &a : vars->CG[1]) {
		double px=a.px;
		double py=a.py;
		double pz=a.pz;
		double v2 = px*px + py*py + pz*pz;
		if (vmax2 < v2) vmax2 = v2 ;
	}
	//double vion2=vars->CG[0][0].px*vars->CG[0][0].px+vars->CG[0][0].py*vars->CG[0][0].py+vars->CG[0][0].pz*vars->CG[0][0].pz;
	//if (vmax2<vion2) vmax2=vion2;
	double vmax = sqrt(vmax2);

	loop_update=margin_length/(vmax*dt)*0.5;
	loop=0;
}

/////////////////////////////////////////////////////////////////////
/*
	make gas-ion pair list
*/
/////////////////////////////////////////////////////////////////////
void
MD::updateInOut(int molFlag){
	for (auto &a : vars->CG[molFlag]){
		double dx = a.qx - vars->CG[0][0].qx;
		double dy = a.qy - vars->CG[0][0].qy;
		double dz = a.qz - vars->CG[0][0].qz;
		double r2 = (dx * dx + dy * dy + dz * dz);
		if (r2 < RI2 && a.inFlag==0) {
			if (molFlag==1) makeDiatomicProp_in(a);
			if (molFlag==2) makePolyatomicProp_in(a);
			a.inFlag=1;
		}
		if (r2 > RI2 && a.inFlag==1) {
			a.inFlag=0;
			if(molFlag==2){
				sprintf(filepath, "vapor_out_%d.dat", int(calculation_number));
				FILE*f=fopen(filepath, "a");
				Molecule *ion=vars->CG[0].data();
				fprintf(f, "%d %e %e %e %e %e %e %e\n", a.id,itime*dt,dx,dy,dz,\
				a.px-vars->CG[0][0].px,a.py-vars->CG[0][0].py,a.pz-vars->CG[0][0].pz);
				fclose(f);
			}
		}
	}
}


/////////////////////////////////////////////////////////////////////
/*
	make vapor-ion pair list
*/
/////////////////////////////////////////////////////////////////////

void
MD::make_pair_gasion(void){
	vars->pairs_gi.clear();
	int i=0;
	for (auto &gin : vars->CG[1]){
		if(gin.inFlag==0) continue;
		int j=0;
		for(auto &ion : vars->CG[0][0].inAtoms){
			double dx=gin.qx-ion.qx;
			double dy=gin.qy-ion.qy;
			double dz=gin.qz-ion.qz;
			double r2 = (dx * dx + dy * dy + dz * dz);
			if(r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairs_gi.push_back(p);
			}
			j++;
		}
		i++;
	}
}

// this may not be needed
void
MD::make_pair_vaporvapor(void){
	vars->pairs_vv.clear();
	const int vs = vars->CG[2].size();
	if(vs>1){
		Pair p;
		for(int i1=0; i1<vs-1; i1++){
			p.i=i1;
			for(int i2=i1+1; i2<vs; i2++){
				p.j=i2;
				vars->pairs_vv.push_back(p);
			}
		}
	}
}

void
MD::make_pair_gasvapor(void){
	vars->pairs_gv.clear();
	int i=0;
	for (auto &gin : vars->CG[1]){
		int j=0;
		for (auto &vin : vars->CG[2]){
			double dx = gin.qx - vin.qx;
			double dy = gin.qy - vin.qy;
			double dz = gin.qz - vin.qz;
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairs_gv.push_back(p);
			}
			j++;
		}
		i++;
	}
}


/////////////////////////////////////////////////////////////////////
/*
	make gas-gas pair list
*/
/////////////////////////////////////////////////////////////////////
void
MD::make_pair_gasgas(void){
	vars->pairs_gg.clear();
	Molecule *gases = vars->CG[1].data();
	int gs=vars->CG[1].size();
	for (int i=0; i<gs-1; i++){
		for (int j=i+1; j<gs; j++){
			double dx = gases[i].qx - gases[j].qx;
			double dy = gases[i].qy - gases[j].qy;
			double dz = gases[i].qz - gases[j].qz;
			adjust_periodic(dx, dy, dz, d_size);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairs_gg.push_back(p);
			}
		}
	}
}


/////////////////////////////////////////////////////////////////////
/*
	check necessity of the pair list updating
*/
/////////////////////////////////////////////////////////////////////
void
MD::check_pairlist(void){
	vars->times.tpair-=omp_get_wtime();
	loop++;
	if(loop>loop_update){
		//boundary_scaling_ion_move();
		boundary_scaling_gas_move();
		boundary_scaling_vapor_move();
		make_pair();
		pre_ion[0]=vars->CG[0][0].qx;
		pre_ion[1]=vars->CG[0][0].qy;
		pre_ion[2]=vars->CG[0][0].qz;
	}
//	if(flags->force_sw==1) sw->check_pairlist(vars);
//	if(flags->force_ters==1) ters->check_pairlist(vars);
	vars->times.tpair+=omp_get_wtime();
}
