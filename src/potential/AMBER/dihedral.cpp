#include "potentialAMBER.hpp"

void
PotentialAMBER::computeDihedral(Variables *vars, FLAG *flags) {
	Atom *ions = vars->Molecules[0].inAtoms.data();
	Dihedral_type *dtypes = vars->dtypes.data();
	Dihedral *dihedrals=vars->Molecules[0].dihedrals.data();
	int dsize=vars->Molecules[0].dihedrals.size();
	for (int idi=0;idi<dsize;idi++) {
		double ff2[3],ff4[3],ff1[3],ff3[3];

		int i=dihedrals[idi].atom1;
		int j=dihedrals[idi].atom2;
		int k=dihedrals[idi].atom3;
		int l=dihedrals[idi].atom4;
		int type=dihedrals[idi].type;

		// 1st bond
		double vb1x = ions[i].qx - ions[j].qx;
		double vb1y = ions[i].qy - ions[j].qy;
		double vb1z = ions[i].qz - ions[j].qz;
		// 2nd bond
		double vb2x = ions[k].qx - ions[j].qx;
		double vb2y = ions[k].qy - ions[j].qy;
		double vb2z = ions[k].qz - ions[j].qz;
		double vb2xm = -vb2x;
		double vb2ym = -vb2y;
		double vb2zm = -vb2z;
		// 3rd bond
		double vb3x = ions[l].qx - ions[k].qx;
		double vb3y = ions[l].qy - ions[k].qy;
		double vb3z = ions[l].qz - ions[k].qz;

		double ax = vb1y*vb2zm - vb1z*vb2ym;
		double ay = vb1z*vb2xm - vb1x*vb2zm;
		double az = vb1x*vb2ym - vb1y*vb2xm;
		double bx = vb3y*vb2zm - vb3z*vb2ym;
		double by = vb3z*vb2xm - vb3x*vb2zm;
		double bz = vb3x*vb2ym - vb3y*vb2xm;
		double rasq = ax*ax + ay*ay + az*az;
		double rbsq = bx*bx + by*by + bz*bz;
		double rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
		double rg = sqrt(rgsq);
		double rginv=0.0;
		double ra2inv=0.0;
		double rb2inv=0.0;
		if (rg > 0) rginv = 1.0/rg;
		if (rasq > 0) ra2inv = 1.0/rasq;
		if (rbsq > 0) rb2inv = 1.0/rbsq;
		double rabinv = sqrt(ra2inv*rb2inv);
		double c = (ax*bx + ay*by + az*bz)*rabinv;
		double s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

		double df = 0.0;

		double p_,df1,ddf1;
		for(int JJ=0;JJ<dtypes[type].multi;JJ++){
			int JJ5=JJ*5;
			p_=1.0;
			ddf1=df1=0.0;
			for (int loop=0; loop<dtypes[type].coeff[JJ5+1]; loop++){
				ddf1 = p_*c - df1*s;
				df1 = p_*s + df1*c;
				p_ = ddf1;
			}
			p_ = p_*dtypes[type].coeff[JJ5+3] + df1*dtypes[type].coeff[JJ5+4];
			df1 = df1*dtypes[type].coeff[JJ5+3] - ddf1*dtypes[type].coeff[JJ5+4];
			df1 *= -dtypes[type].coeff[JJ5+1];
			p_ += 1.0;
	        if(dtypes[type].coeff[1]==0){
	            p_=1.0+dtypes[type].coeff[JJ5+3];
	            df1=0.0;
	        }
			df += (-dtypes[type].coeff[JJ5] * df1);
			if (flags->eflag) vars->Utotal.Uion+= dtypes[type].coeff[JJ5] * p_;
		}


		double fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
		double hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
		double fga = fg*ra2inv*rginv;
		double hgb = hg*rb2inv*rginv;
		double gaa = -ra2inv*rg;
		double gbb = rb2inv*rg;
		double dtfx = gaa*ax;
		double dtfy = gaa*ay;
		double dtfz = gaa*az;
		double dtgx = fga*ax - hgb*bx;
		double dtgy = fga*ay - hgb*by;
		double dtgz = fga*az - hgb*bz;
		double dthx = gbb*bx;
		double dthy = gbb*by;
		double dthz = gbb*bz;
		double sx2 = df*dtgx;
		double sy2 = df*dtgy;
		double sz2 = df*dtgz;
		ff1[0] = df*dtfx;
		ff1[1] = df*dtfy;
		ff1[2] = df*dtfz;
		ff2[0] = sx2 - ff1[0];
		ff2[1] = sy2 - ff1[1];
		ff2[2] = sz2 - ff1[2];
		ff4[0] = df*dthx;
		ff4[1] = df*dthy;
		ff4[2] = df*dthz;
		ff3[0] = -sx2 - ff4[0];
		ff3[1] = -sy2 - ff4[1];
		ff3[2] = -sz2 - ff4[2];

		ions[i].fx += ff1[0];
		ions[i].fy += ff1[1];
		ions[i].fz += ff1[2];
		ions[j].fx += ff2[0];
		ions[j].fy += ff2[1];
		ions[j].fz += ff2[2];
		ions[k].fx += ff3[0];
		ions[k].fy += ff3[1];
		ions[k].fz += ff3[2];
		ions[l].fx += ff4[0];
		ions[l].fy += ff4[1];
		ions[l].fz += ff4[2];

	}
}
