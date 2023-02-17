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

//options
	int fix_cell_center;
	int gyration;
	int RDF;

	bool inter;

//couase-grained
	int cg;

	int eflag;

	FLAG(void);
	~FLAG(void);

private:

};
