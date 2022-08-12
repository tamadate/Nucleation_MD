//------------------------------------------------------------------------
#include "CalculationCond.hpp"
#include "PhysicalProp.hpp"
#include "md.hpp"
//------------------------------------------------------------------------



int main ( int argc,char *argv[] ) {
        
/////////////////////////////////////////////////////////////////////
/*	
	Diffusion coefficient estimation	
	-Diffusion coefficient and collision rate coefficient simulation 
	require 2 and 7 input parameters, respectively.
	- 1:molecular infomation file name 
	- 2:calculation condition file name 
	- 3:calculation number
*/
/////////////////////////////////////////////////////////////////////

	if ( argc==3 ) {
		calculation_number = stoi ( argv[ 2 ] );
		MD *md = new MD ( argv[ 1 ] );	
		md -> mbdist -> makeWeightedMB(md->pp->cgas,md->pp->mgas);
		md -> mbdistV -> makeWeightedMB(md->pp->cvapor,md->pp->mvapor);
		md -> run_diff ( argv );
	}



/*Exception error*/
	if ( argc != 3 ) cout<<"Error:Number of input parameters. -> Diffusion coefficient simulation require 2 input parameters \n1:molecular infomation file name \n2:calculation condition file name \n3:Calculation number"<<endl;
	return 0;
}

