#include "md.hpp"


/*########################################################################################

-----Calculating center of mass and velocity of center of mass-----

#######################################################################################*/

void
MD::getIonCenterProp(void) {
	for ( int i = 0; i < 3; i++ ) vars->IonX[ i ] = 0, vars->IonV[ i ] = 0;
	for ( auto &a : vars -> ions)	{
		vars->IonX[ 0 ] += a.qx * a.mass;
		vars->IonX[ 1 ] += a.qy * a.mass;
		vars->IonX[ 2 ] += a.qz * a.mass;
		vars->IonV[ 0 ] += a.px * a.mass;
		vars->IonV[ 1 ] += a.py * a.mass;
		vars->IonV[ 2 ] += a.pz * a.mass;
	}
	for ( int i = 0; i < 3; i++ ) vars->IonX[ i ] /= pp -> Mion;
	for ( int i = 0; i < 3; i++ ) vars->IonV[ i ] /= pp -> Mion;
}

void
MD::getGasCenterProp(void) {
	for ( int i = 0; i < 3; i++) vars->gasX[ i ] = 0, vars->gasV[ i ] = 0;
	for (auto &a : vars -> gases)	{
		vars->gasX[ 0 ] += a.qx * a.mass;
		vars->gasX[ 1 ] += a.qy * a.mass;
		vars->gasX[ 2 ] += a.qz * a.mass;
		vars->gasV[ 0 ] += a.px *a.mass;
		vars->gasV[ 1 ] += a.py *a.mass;
		vars->gasV[ 2 ] += a.pz *a.mass;
	}
	for ( int i = 0; i < 3; i++ ) vars->gasX[ i ] /= pp -> Mgas;
	for ( int i = 0; i < 3; i++) vars->gasV[ i ] /= pp -> Mgas;
}


