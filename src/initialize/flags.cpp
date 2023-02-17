//------------------------------------------------------------------------
#include "../flags.hpp"

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
	cg=0;
	eflag=0;
	inter=false;
}

/**********************************destractor******************************************/
FLAG::~FLAG(void) {

}
