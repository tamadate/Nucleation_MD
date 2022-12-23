#include "../md.hpp"

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
