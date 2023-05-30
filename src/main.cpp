#include "constants.hpp"
#include "PhysicalProp.hpp"
#include "md.hpp"

/////////////////////////////////////////////////////////////////////
/*
	Title:
		Vapor uptake molecular dynamics (MD) simulation
	Author:
		Dr. Tomoya Tamadate 
	Affiliation:
		University of Minnesota
		Kanazawa University
	Code development was started in 2021 with Chris Hogan (UMN)
*/
/////////////////////////////////////////////////////////////////////

int main ( int argc,char *argv[] ) {

	// If two input parameters are feeded properly, the calculation is stopped
	// input 1:molecular infomation file name
	// input 2:starting calculation number
	if(argc==3){
		// setting OpenMP for parallel computing
		int Nth=omp_get_max_threads(); 
		omp_set_num_threads(Nth);

		// main part
		#pragma omp parallel for
		for(int i=0;i<Nth;i++){
			// generate MD class
			MD *md=new MD(argv[1],stoi(argv[2])+i);
			// run MD simulation
			md->run(argv);
		}
	}

	// in case of error, it dump an error message
	if (argc!=3) cout<<"Error: Number of input parameters.\n \
		-> Diffusion coefficient simulation require 2 input parameters \n \
		1: Calculation condition file name \n \
		2: Starting calculation number"<<endl;

	return 0;
}
