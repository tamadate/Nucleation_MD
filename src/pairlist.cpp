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
	make_pair_gasion();
	make_pair_vaporion();
	make_pair_gasvapor();
	if(flags->inter_gg==1) make_pair_gasgas();
//	set number of steps to update pair list
	Atom *gases = vars->gases.data();
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
MD::make_pair_gasion(void){
// clear vars->gas_in, vars->gas_out and pair list of gas-ion
	vars->gas_in.clear();
	vars->gas_out.clear();
	vars->pairs_gi.clear();

	Atom *gases = vars->gases.data();
	Atom *ions = vars->ions.data();
	int is=vars->ions.size();
	for (int i=0;i<Nof_around_gas;i++){
		int flag_in=0;
		for (int j=0;j<is;j++){
			double dx = gases[i].qx - ions[j].qx;
			double dy = gases[i].qy - ions[j].qy;
			double dz = gases[i].qz - ions[j].qz;
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				flag_in+=1;
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairs_gi.push_back(p);
                if(pp->num_gas==2) {
                    p.i=i+Nof_around_gas;
                    vars->pairs_gi.push_back(p);
                }
			}
		}
//if inter-gas interaction flag is ON, stand flag_in ON
//if flag_in ON, molecule push back to vars->gas_in
//if flag_in OFF, molecule push back to vars->gas_out
		if(flags->inter_gg==1) flag_in=1;	
		if(flag_in>0) {
            vars->gas_in.push_back(i);
            if(pp->num_gas==2) {
               makeDiatomicProp_in(i);
               vars->gas_in.push_back(i+Nof_around_gas);
            }
        }
		if(flag_in==0) {
            vars->gas_out.push_back(i);
            if(pp->num_gas==2) {
                makeDiatomicProp_out(i);
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
MD::make_pair_vaporion(void){
// clear vars->gas_in, vars->gas_out and pair list of gas-ion
	Molecule *vapors = vars->vapors.data();
	for (auto i : vars->vapor_in){
		double X=0;
		double Y=0;
		double Z=0;
		double VX=0;
		double VY=0;
		double VZ=0;
		double Mass=0;
		for (auto &a : vars->vapors[i].inAtoms){
			X+=a.qx*a.mass;
			Y+=a.qy*a.mass;
			Z+=a.qz*a.mass;
			VX+=a.px*a.mass;
			VY+=a.py*a.mass;
			VZ+=a.pz*a.mass;
			Mass+=a.mass;
		}
        //  averaged position
	    vars->vapors[i].qx=X/Mass;
	    vars->vapors[i].qy=Y/Mass;
	    vars->vapors[i].qz=Z/Mass;
        //  averaged velocity
        vars->vapors[i].px=VX/Mass;
        vars->vapors[i].py=VY/Mass;
        vars->vapors[i].pz=VZ/Mass;
	}

	vars->vapor_in.clear();
	vars->vapor_out.clear();

	for (int i=0;i<Nof_around_vapor;i++){
		int flag_in=0;
		for (auto &a : vars->ions){
			double dx = vapors[i].qx - a.qx;
			double dy = vapors[i].qy - a.qy;
			double dz = vapors[i].qz - a.qz;
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				flag_in=1;
				break;
			}
		}
		if(flag_in>0) {
			vars->vapor_in.push_back(i);
			makePolyatomicProp_in(i);
		}
		if(flag_in==0) {
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
MD::make_pair_gasvapor(void){
// clear vars->gas_in, vars->gas_out and pair list of gas-ion
	Molecule *vapors = vars->vapors.data();
	Atom *gases = vars->gases.data();
	Pair p;
	for (auto i : vars->gas_in){
		for (auto j : vars->vapor_in){
			double dx = gases[i].qx - vapors[j].qx;
			double dy = gases[i].qy - vapors[j].qy;
			double dz = gases[i].qz - vapors[j].qz;
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				p.i=i;
				p.j=j;
				pairs_gv.push_back(p);
                if(pp->num_gas==2) {
                    p.i=i+Nof_around_gas;
                    pairs_gv.push_back(p);
                }
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
	Atom *gases = vars->gases.data();
	int gs=vars->gases.size();
	for (int i=0; i<gs-1; i++){
		for (int j=i+1; j<gs; j++){
			double dx = gases[i].qx - gases[j].qx;
			double dy = gases[i].qy - gases[j].qy;
			double dz = gases[i].qz - gases[j].qz;
			adjust_periodic(dx, dy, dz);
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
	loop++;
	if(loop>loop_update){
		Atom *gases = vars->gases.data();
		Molecule *vapors = vars->vapors.data();
		boundary_scaling_ion_move();
		for (auto &i : vars->gas_out) {
            gases[i].qx += gases[i].px*dt*loop;
            gases[i].qy += gases[i].py*dt*loop;
            gases[i].qz += gases[i].pz*dt*loop;
		}
		for (auto &i : vars->vapor_out) {
            vapors[i].qx += vapors[i].px*dt*loop;
            vapors[i].qy += vapors[i].py*dt*loop;
            vapors[i].qz += vapors[i].pz*dt*loop;
		}
		boundary_scaling_gas_move();
		boundary_scaling_vapor_move();
		make_pair();
		for(int i=0;i<3;i++) pre_ion[i]=ion_r[i];
	}
//	if(flags->force_sw==1) sw->check_pairlist(vars);
//	if(flags->force_ters==1) ters->check_pairlist(vars);
}




