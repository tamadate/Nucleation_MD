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


	FLAG(void);
	~FLAG(void);

private:

};
