#include "../md.hpp"


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
	update_gas_in();
	update_vapor_in();
	for (auto &a : InterInter) a->makePair(vars);
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
	double vion2=vars->IonV[0]*vars->IonV[0]+vars->IonV[1]*vars->IonV[1]+vars->IonV[2]*vars->IonV[2];
	if (vmax2<vion2) vmax2=vion2;
	double vmax = sqrt(vmax2);

	loop_update=con->MARGIN/(vmax*dt);
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
// clear vars->gas_in, vars->gas_out and pair list of gas-ion
	vars->gas_in.clear();
	vars->gas_out.clear();

	for (int i=0;i<con->Nof_around_gas;i++){
		double dx,dy,dz;
		double r2 = distFromIonCenter(gases[i],dx,dy,dz);
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

	vars->vapor_in.clear();
	vars->vapor_out.clear();

	for (int i=0;i<con->Nof_around_vapor;i++){
		double dx,dy,dz;
		double r2 = distFromIonCenter(vapors[i],dx,dy,dz);
		if (r2 < RI2) {
			if(r2<100 && collisionFlagVapor[i]==0){
				obs->Ovin(i,itime*dt);
				collisionFlagVapor[i]=itime;
			}
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
		double coeff=dt*loop;
		for (auto &i : vars->gas_out) {
			gases[i].qx += gases[i].px*coeff;
			gases[i].qy += gases[i].py*coeff;
			gases[i].qz += gases[i].pz*coeff;
		}
		for (auto &i : vars->vapor_out) {
			vapors[i].qx += vapors[i].px*coeff;
			vapors[i].qy += vapors[i].py*coeff;
			vapors[i].qz += vapors[i].pz*coeff;
		}
		boundary_scaling_gas_move();
		boundary_scaling_vapor_move();
		boundary_scaling_ion_move();
		make_pair();
		vars->preIonX[0]=vars->IonX[0];
		vars->preIonX[1]=vars->IonX[1];
		vars->preIonX[2]=vars->IonX[2];
	}
//	if(flags->force_sw==1) sw->check_pairlist(vars);
//	if(flags->force_ters==1) ters->check_pairlist(vars);
	vars->times.tpair+=omp_get_wtime();
}
