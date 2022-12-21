# pragma once

class FLAG{
public:
//NVT flags
	int velocity_scaling;
	int nose_hoover_ion;
	int nose_hoover_gas;
	int semi_NVT_gasgas;

//dump flags
	int dump_gas;
	int dump_fix;
	int dump_2nd;

//intra-atomic interaction in ion
	int force_ters;
	int force_sw;
	int intra_AMBER;
	int force_born;

//inter-atomic ion-gas
	int force_lj;
	int force_ion_dipole;

//inter-atomic gas-gas
	int inter_vi;

//inter-atomic gas-gas
	int inter_gg;

//inter-atomic vapor-vapor
	int inter_vv;

//inter-atomic vapor-gas
	int inter_vg;

//intra-atomic vapor
	int vapor_intra;

//intra-atomic gas
	int gas_intra;

//Electric field
	int efield;

//options
	int fix_cell_center;
	int gyration;
	int RDF;

//couase-grained
	int cg;

	int eflag;
	int takeOver;


	FLAG(void){
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
		takeOver=0;

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

		inter_vv=0;

		inter_vg=0;

		inter_vi=0;

		vapor_intra=1;

		gas_intra=1;

	//couase-grain
		cg=0;

		eflag=0;
	}

	~FLAG(void){};

private:

};
