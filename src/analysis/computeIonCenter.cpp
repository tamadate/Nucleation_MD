//------------------------------------------------------------------------
#include "../md.hpp"

void
MD::analysis_ion(void) {
	Molecule *ion = vars -> Molecules.data();
	ion[0].qx=0;
	ion[0].qy=0;
	ion[0].qz=0;
	ion[0].px=0;
	ion[0].py=0;
	ion[0].pz=0;

	for ( auto &a : ion[0].inAtoms)	{
		ion[0].qx += a.qx * a.mass;
		ion[0].qy += a.qy * a.mass;
		ion[0].qz += a.qz * a.mass;
		ion[0].px += a.px * a.mass;
		ion[0].py += a.py * a.mass;
		ion[0].pz += a.pz * a.mass;
	}
	ion[0].qx /= pp -> Mion;
	ion[0].qy /= pp -> Mion;
	ion[0].qz /= pp -> Mion;
	ion[0].px /= pp -> Mion;
	ion[0].py /= pp -> Mion;
	ion[0].pz /= pp -> Mion;
}
