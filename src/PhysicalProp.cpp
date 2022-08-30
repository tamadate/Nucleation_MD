//------------------------------------------------------------------------
#include "PhysicalProp.hpp"
//------------------------------------------------------------------------


#define SIZE 1000000

PhysicalProp::PhysicalProp(void){
	gastype=2;	/*1:He, 2:Ar, 3:N2*/
	vaportype=1;	/*1:He, 2:Ar, 3:N2*/
	Noftimestep=10000000000;
	p=pow(10,(2.0+3.0));	/*Pa*/
	T = 300;     /* temp. */
	Nof_around_gas=10000;
	Nof_around_vapor=1;
}
