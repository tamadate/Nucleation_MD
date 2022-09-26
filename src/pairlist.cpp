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
	update_gas_in();
	update_vapor_in();
	make_pair_gasion();
	make_pair_gasvapor();
	make_pair_vaporvapor();
	if(flags->inter_gg==1) make_pair_gasgas();
//	set number of steps to update pair list
	Molecule *gases = vars->gases.data();
	double vmax2 = 0.0;
	int gos=vars->gas_out.size();
	for (auto &a : vars->gas_out) {
		double px=gases[a].px;
		double py=gases[a].py;
		double pz=gases[a].pz;
		double v2 = px*px + py*py + pz*pz;
		if (vmax2 < v2) vmax2 = v2;
	}
	double vion2=ion_v[0]*ion_v[0]+ion_v[1]*ion_v[1]+ion_v[2]*ion_v[2];
	if (vmax2<vion2) vmax2=vion2;
	double vmax = sqrt(vmax2);

	loop_update=margin_length/(vmax*dt);
	loop=0;
}

/////////////////////////////////////////////////////////////////////
/*
	make gas-ion pair list
*/
/////////////////////////////////////////////////////////////////////
void
MD::update_gas_in(void){
	Molecule *gases = vars->gases.data();
	updateGasinCenters();
// clear vars->gas_in, vars->gas_out and pair list of gas-ion
	vars->gas_in.clear();
	vars->gas_out.clear();


	for (int i=0;i<Nof_around_gas;i++){
		double dx = gases[i].qx - ion_r[0];
		double dy = gases[i].qy - ion_r[1];
		double dz = gases[i].qz - ion_r[2];
		double r2 = (dx * dx + dy * dy + dz * dz);
		if (r2 < RI2){
//if inter-gas interaction flag is ON, stand flag_in ON
//if flag_in ON, molecule push back to vars->gas_in
//if flag_in OFF, molecule push back to vars->gas_out
    	vars->gas_in.push_back(i);
    	makeDiatomicProp_in(i);
	  }
		else{
      vars->gas_out.push_back(i);
      makeDiatomicProp_out(i);
		}
	}
}

/////////////////////////////////////////////////////////////////////
/*
	make vapor-ion pair list
*/
/////////////////////////////////////////////////////////////////////
void
MD::update_vapor_in(void){
	Molecule *vapors = vars->vapors.data();
// clear vars->gas_in, vars->gas_out and pair list of gas-ion
	updateVaporinCenters();

	vars->vapor_in.clear();
	vars->vapor_out.clear();

	for (int i=0;i<Nof_around_vapor;i++){
		double dx = vapors[i].qx - ion_r[0];
		double dy = vapors[i].qy - ion_r[1];
		double dz = vapors[i].qz - ion_r[2];
		double r2 = (dx * dx + dy * dy + dz * dz);
		if (r2 < RI2) {
			vars->vapor_in.push_back(i);
			makePolyatomicProp_in(i);
		}
		else{
			vars->vapor_out.push_back(i);
			makePolyatomicProp_out(i);
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
	Molecule *gases = vars->gases.data();
	Atom *ions = vars->ions.data();
	for (auto i : vars->gas_in){
		for(int j=0;j<vars->ions.size();j++){
			double dx=gases[i].qx-ions[j].qx;
			double dy=gases[i].qy-ions[j].qy;
			double dz=gases[i].qz-ions[j].qz;
			double r2 = (dx * dx + dy * dy + dz * dz);
			if(r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairs_gi.push_back(p);
			}
		}
	}
}

void
MD::make_pair_vaporvapor(void){
	vars->pairs_vv.clear();
	Molecule *vapors = vars->vapors.data();
	const int vs = vars->vapor_in.size();
	if(vs>1){
		Pair p;
		for(int i1=0; i1<vs-1; i1++){
			p.i=vars->vapor_in[i1];
			for(int i2=i1+1; i2<vs; i2++){
				p.j=vars->vapor_in[i2];
				vars->pairs_vv.push_back(p);
			}
		}
	}
}

void
MD::make_pair_gasvapor(void){
	vars->pairs_gv.clear();
	Molecule *vapors = vars->vapors.data();
	Molecule *gases = vars->gases.data();
	for (auto i : vars->gas_in){
		for (auto j : vars->vapor_in){
			double dx = gases[i].qx - vapors[j].qx;
			double dy = gases[i].qy - vapors[j].qy;
			double dz = gases[i].qz - vapors[j].qz;
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairs_gv.push_back(p);
			}
		}
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
	Molecule *gases = vars->gases.data();
	int gs=vars->gases.size();
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
		Molecule *gases = vars->gases.data();
		Molecule *vapors = vars->vapors.data();
		boundary_scaling_ion_move();
		//#pragma omp parallel for
		for (auto &i : vars->gas_out) {
        gases[i].qx += gases[i].px*dt*loop;
        gases[i].qy += gases[i].py*dt*loop;
        gases[i].qz += gases[i].pz*dt*loop;
		}
		//#pragma omp parallel for
		for (auto &i : vars->vapor_out) {
        vapors[i].qx += vapors[i].px*dt*loop;
        vapors[i].qy += vapors[i].py*dt*loop;
        vapors[i].qz += vapors[i].pz*dt*loop;
		}
		boundary_scaling_gas_move();
		boundary_scaling_vapor_move();
		make_pair();
		pre_ion[0]=ion_r[0];
		pre_ion[1]=ion_r[1];
		pre_ion[2]=ion_r[2];
	}
//	if(flags->force_sw==1) sw->check_pairlist(vars);
//	if(flags->force_ters==1) ters->check_pairlist(vars);
	vars->times.tpair+=omp_get_wtime();
}
