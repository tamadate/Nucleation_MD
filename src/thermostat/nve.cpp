//------------------------------------------------------------------------
#include "../md.hpp"
//------------------------------------------------------------------------

void
MD::setNVE(void){
	flags->velocity_scaling=0;
	flags->nose_hoover_ion=0;
	flags->nose_hoover_gas=0;
}
