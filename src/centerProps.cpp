#include "../md.hpp"


/*########################################################################################

-----Calculating center of mass and velocity of center of mass-----

#######################################################################################*/

void
MD::getIonCenterProp(void) {
	for ( int i = 0; i < 3; i++ ) vars->ion_r[ i ] = 0, vars->ion_v[ i ] = 0;
	for ( auto &a : vars -> ions)	{
		vars->ion_r[ 0 ] += a.qx * a.mass;
		vars->ion_r[ 1 ] += a.qy * a.mass;
		vars->ion_r[ 2 ] += a.qz * a.mass;
		vars->ion_v[ 0 ] += a.px * a.mass;
		vars->ion_v[ 1 ] += a.py * a.mass;
		vars->ion_v[ 2 ] += a.pz * a.mass;
	}
	for ( int i = 0; i < 3; i++ ) vars->ion_r[ i ] /= pp -> Mion;
	for ( int i = 0; i < 3; i++ ) vars->ion_v[ i ] /= pp -> Mion;
}

void
MD::getGasCenterProp(void) {
	for ( int i = 0; i < 3; i++) vars->gas_r[ i ] = 0, vars->gas_v[ i ] = 0;
	for (auto &a : vars -> gases)	{
		vars->gas_r[ 0 ] += a.qx * a.mass;
		vars->gas_r[ 1 ] += a.qy * a.mass;
		vars->gas_r[ 2 ] += a.qz * a.mass;
		vars->gas_v[ 0 ] += a.px *a.mass;
		vars->gas_v[ 1 ] += a.py *a.mass;
		vars->gas_v[ 2 ] += a.pz *a.mass;
	}
	for ( int i = 0; i < 3; i++ ) vars->gas_r[ i ] /= pp -> Mgas;
	for ( int i = 0; i < 3; i++) vars->gas_v[ i ] /= pp -> Mgas;
}


