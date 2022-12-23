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
	updateInOut();
	make_pair_gasion();
	make_pair_short();
	make_pair_vaporvapor();
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

//------------------------------------------------------------//
/*	make gas-ion pair list  */
//------------------------------------------------------------//
void
MD::updateInOut(void){
	Molecule *mols=vars->Molecules.data();
	for (auto i : vars->MolID[1]){
		double dx = mols[i].qx - mols[0].qx;
		double dy = mols[i].qy - mols[0].qy;
		double dz = mols[i].qz - mols[0].qz;
		double r2 = (dx * dx + dy * dy + dz * dz);
		double dr2=dx*dx+dy*dy+dz*dz;
		int oriFlag=Region[i];
		if(dr2<RO2) Region[i]=AA;
		else {
			if(dr2<RI2) {
				Region[i]=AACG;
				if(oriFlag==CG) makeDiatomicProp_in(mols[i]);
			}
			else Region[i]=CG;
		}
	}
	for (auto i : vars->MolID[2]){
		double dx = mols[i].qx - mols[0].qx;
		double dy = mols[i].qy - mols[0].qy;
		double dz = mols[i].qz - mols[0].qz;
		double r2 = (dx * dx + dy * dy + dz * dz);
		double dr2=dx*dx+dy*dy+dz*dz;
		int oriFlag=Region[i];
		if(dr2<RO2) Region[i]=AA;
		else {
			if(dr2<RI2) {
				Region[i]=AACG;
				if(oriFlag==CG) makePolyatomicProp_in(mols[i]);
			}
			else {
				Region[i]=CG;
				if(oriFlag!=CG){
					sprintf(filepath, "vapor_out_%d.dat", int(calculation_number));
					FILE*f=fopen(filepath, "a");
					fprintf(f, "%d %e %e %e %e %e %e %e\n", mols[0].id,itime*dt,dx,dy,dz,\
					mols[i].px-mols[0].px,mols[i].py-mols[0].py,mols[i].pz-mols[0].pz);
					fclose(f);
				}
			}
		}
	}
}


//------------------------------------------------------------//
/*	Make ion-gas interaction pair list for gas molecule  */
//------------------------------------------------------------//

void
MD::make_pair_gasion(void){
	vars->pairs_gi.clear();
	Molecule *mols=vars->Molecules.data();
	for (auto &i : vars->MolID[1]){
		if(vars->Region[i]==CG) continue;
		int j=0;
		for(auto &ion : vars->Molecules[0].inAtoms){
			double dx=mols[i].qx-ion.qx;
			double dy=mols[i].qy-ion.qy;
			double dz=mols[i].qz-ion.qz;
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


//------------------------------------------------------------//
/*	Make short range interaction pair list for gas molecule  */
//------------------------------------------------------------//
void
MD::make_pair_short(void){
	vars->pairsShort.clear();
	Molecule *mols=vars->Molecules.data();

	// make ion-gas pair list
	for (auto i : vars->MolID[1]){
		for (auto j : vars->MolID[2]){
			double dx = mols[i].qx - mols[j].qx;
			double dy = mols[i].qy - mols[j].qy;
			double dz = mols[i].qz - mols[j].qz;
			adjust_periodic(dx, dy, dz, d_size);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairsShort.push_back(p);
			}
		}
	}

	// make gas-gas pair list
	int gs=vars->MolID[1].size();
	for (int i=0; i<gs-1; i++){
		for (int j=i+1; j<gs; j++){
			double dx = mols[i].qx - mols[j].qx;
			double dy = mols[i].qy - mols[j].qy;
			double dz = mols[i].qz - mols[j].qz;
			adjust_periodic(dx, dy, dz, d_size);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairsShort.push_back(p);
			}
		}
	}
}

//------------------------------------------------------------//
/*	Make vapor-vapor interaction pair list for gas molecule  */
//------------------------------------------------------------//
void
MD::make_pair_vaporvapor(void){
	vars->pairs_vv.clear();
	Molecule *mols=vars->Molecules.data();

	// make gas-gas pair list
	int vs=vars->MolID[2].size();
	for (int i=0; i<vs-1; i++){
		for (int j=i+1; j<vs; j++){
			double dx = mols[i].qx - mols[j].qx;
			double dy = mols[i].qy - mols[j].qy;
			double dz = mols[i].qz - mols[j].qz;
			adjust_periodic(dx, dy, dz, d_size);
			double r2 = (dx * dx + dy * dy + dz * dz);
			// 00000001|00000001=00000001
			// 00000011|00000001=00000011
			// 00000011|00000011=00000011
			// 00000010|00000011=00000011
			// 00000010|00000010=00000010 ... only this case, no AA interaction

			if((Region[i]|Region[j])==CG){	// if both of two are CG region
				if (r2 < ML2){
					Pair p;
					p.i=i;
					p.j=j;
					vars->pairsShort.push_back(p);
				}
			}
			else{
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairs_vv.push_back(p);
			}
		}
	}
}


//------------------------------------------------------------//
/*	check necessity of the pair list updating	*/
//------------------------------------------------------------//
void
MD::check_pairlist(void){
	vars->times.tpair-=omp_get_wtime();
	loop++;
	if(loop>loop_update){
		make_pair();
		pre_ion[0]=vars->Molecules[0].qx;
		pre_ion[1]=vars->Molecules[0].qy;
		pre_ion[2]=vars->Molecules[0].qz;
	}
//	if(flags->force_sw==1) sw->check_pairlist(vars);
//	if(flags->force_ters==1) ters->check_pairlist(vars);
	vars->times.tpair+=omp_get_wtime();
}
