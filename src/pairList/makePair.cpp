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
	updateInOut();
	make_pairLJ();
  make_pairLJCoul();
	make_pairLJHybrid();
	make_pairsLJCoulHybrid();
	// set number of steps for next pair list update
	double vmax2 = 0.0;
	for (auto i : vars->MolID[1]) {
		double px=vars->Molecules[i].px;
		double py=vars->Molecules[i].py;
		double pz=vars->Molecules[i].pz;
		double v2 = px*px + py*py + pz*pz;
		if (vmax2 < v2) vmax2 = v2 ;
	}
	//double vion2=vars->Molecules[0][0].px*vars->Molecules[0][0].px+vars->Molecules[0][0].py*vars->Molecules[0][0].py+vars->Molecules[0][0].pz*vars->Molecules[0][0].pz;
	//if (vmax2<vion2) vmax2=vion2;
	double vmax = sqrt(vmax2);

	loop_update=margin_length/(vmax*dt)*0.5;
	loop=0;
}
