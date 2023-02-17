//------------------------------------------------------------------------
#include "constants.hpp"
#include "PhysicalProp.hpp"
#include "md.hpp"
//------------------------------------------------------------------------



int main ( int argc,char *argv[] ) {

/////////////////////////////////////////////////////////////////////
/*
	Single ion tracking
	-Diffusion coefficient and vapor uptaking simulations
	require 2 input parameters.
	- 1:molecular infomation file name
	- 2:calculation number
*/
/////////////////////////////////////////////////////////////////////

	if(argc==3){
		int Nth=omp_get_max_threads();
		omp_set_num_threads(Nth);
		#pragma omp parallel for
		for(int i=0;i<Nth;i++){
			MD *md=new MD(argv[1],stoi(argv[2])+i);
			md->run(argv);
		}
	}
	if (argc!=3) cout<<"Error:Number of input parameters. -> Diffusion coefficient simulation require 2 input parameters \n1:molecular infomation file name \n2:calculation condition file name \n3:Calculation number"<<endl;
	return 0;
}
