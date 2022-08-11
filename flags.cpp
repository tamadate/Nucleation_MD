//------------------------------------------------------------------------
#include "md.hpp"

/**********************************constructor******************************************/
FLAG::FLAG(void) {
	velocity_scaling=0;
	nose_hoover_ion=0;
	nose_hoover_gas=0;
	semi_NVT_gasgas=0;
	dump_gas=1;
	dump_fix=1;
	fix_cell_center=0;
	dump_2nd=1;
	gyration=0;
	RDF=0;

//	ion intratomic interaction
	force_ters=0;
	force_sw=0;
	force_born=0;
	intra_AMBER=0;

//	ion-gas interaction
	force_lj=0;
	force_ion_dipole=0;

//	Electric field
	efield=0;

//	gas-gas interaction
	inter_gg=0;

	inter_vv=1;

	inter_vg=1;

	inter_vi=1;

	vapor_intra=1;

	gas_intra=1;

//couase-grain
	cg=0;

	eflag=0;


}

/**********************************destractor******************************************/
FLAG::~FLAG(void) {

}


