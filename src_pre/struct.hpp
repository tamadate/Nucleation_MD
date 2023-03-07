
struct Atom_type {
	double mass;
	std::string name;
	double coeff1, coeff2;
};
//------------------------------------------------------------------------
struct Bond_type {
	double coeff[2];
};
//------------------------------------------------------------------------
struct Angle_type {
	double coeff[2];
};
//------------------------------------------------------------------------
struct Dihedral_type {
 	int multi;
	double coeff[15];
};
//------------------------------------------------------------------------
struct Atom {
 	double qx, qy, qz;
	double px, py, pz;
	double charge;
	int type;
	int id;
	double fx, fy, fz;
	double mass;
	int ix, iy, iz;
};

//------------------------------------------------------------------------
struct Bond {
	int atom1, atom2, type;
};
//------------------------------------------------------------------------
struct Angle {
	int atom1, atom2, atom3, type;
};
//------------------------------------------------------------------------
struct Dihedral {
	int atom1, atom2, atom3, atom4, type;
};
//------------------------------------------------------------------------
struct Pair{
	int i,j;
};
struct Pair_many{
	int i;
	std::vector<int> j;	/*	pair list	*/
};

//------------------------------------------------------------------------
struct N2_relative{
	double x,y,z,vx1,vy1,vz1,vx2,vy2,vz2;
};
//------------------------------------------------------------------------
struct Molecule {
 	double qx, qy, qz;
	double px, py, pz;
	double mass;
	int ix, iy, iz;
	std::vector<Atom> inAtoms;
	std::vector<Bond> bonds;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;
	int inFlag;
};

struct Potentials {
	double Uion;
	double Ugas;
	double Uvap;
	double Ugi;
	double Ugg;
	double Uvg;
	double Uvi;
	double Uvv;
};

struct Times {
  double tion;
  double tgas;
  double tgi;
  double tvap;
  double tvv;
  double tvg;
  double tvi;
  double tpair;
  double tvel;
  double tpos;
  double tetc;
};
