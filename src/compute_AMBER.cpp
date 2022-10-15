#include "potential.hpp"
/*########################################################################################

-----compute intramolecular interaction-----

#######################################################################################*/

/**********************************Force calculation******************************************/
void
PotentialAMBER::compute(Variables *vars, FLAG *flags) {
	vars->times.tion-=omp_get_wtime();
	computeLong(vars,flags);
	computeBond(vars,flags);
	computeAngle(vars,flags);
	computeDihedral(vars,flags);
	vars->times.tion+=omp_get_wtime();
}


void
PotentialAMBER::computeLong(Variables *vars, FLAG *flags) {
	Atom *ions = vars->effectiveIn[0][0].inAtoms.data();
	int lpsize=longPair.size();
	#pragma omp parallel for
	for (int ip=0;ip<lpsize;ip++) {
		int nth=omp_get_thread_num();
		int i=longPair[ip].i;
		int j=longPair[ip].j;
		double dx = ions[i].qx - ions[j].qx;
		double dy = ions[i].qy - ions[j].qy;
		double dz = ions[i].qz - ions[j].qz;
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r2inv = 1/rsq;
		int type1=ions[i].type;
		int type2=ions[j].type;
		double r6inv = r2inv * r2inv * r2inv;
		double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
		double force_coul = qqrd2e * ions[i].charge * ions[j].charge * sqrt(r2inv);
		double force_pair = (force_lj + force_coul)*r2inv;
		ions[i].fxMP[nth] += force_pair * dx;
		ions[i].fyMP[nth] += force_pair * dy;
		ions[i].fzMP[nth] += force_pair * dz;
		ions[j].fxMP[nth] -= force_pair * dx;
		ions[j].fyMP[nth] -= force_pair * dy;
		ions[j].fzMP[nth] -= force_pair * dz;
		if(flags->eflag) {
			vars->U_MP[nth].Uion+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			vars->U_MP[nth].Uion+=force_coul;
		}
	}
}

void
PotentialAMBER::computeBond(Variables *vars, FLAG *flags) {
	Atom *ions = vars->effectiveIn[0][0].inAtoms.data();
	Bond_type *btypes = vars->btypes.data();
	Bond *bonds=vars->effectiveIn[0][0].bonds.data();
	int bsize=vars->effectiveIn[0][0].bonds.size();
	#pragma omp parallel for
	for (int ib=0;ib<bsize;ib++) {
		int nth=omp_get_thread_num();
		int i=bonds[ib].atom1;
		int j=bonds[ib].atom2;
		int type=bonds[ib].type;
		double dx = ions[i].qx - ions[j].qx;
		double dy = ions[i].qy - ions[j].qy;
		double dz = ions[i].qz - ions[j].qz;
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r = sqrt(rsq);
		double dr = (r-btypes[type].coeff[1]);
		double rk = btypes[type].coeff[0] * dr;
		double force_bond_harmonic;
		force_bond_harmonic = -2.0*rk/r;
		ions[i].fxMP[nth] += force_bond_harmonic * dx;
		ions[i].fyMP[nth] += force_bond_harmonic * dy;
		ions[i].fzMP[nth] += force_bond_harmonic * dz;
		ions[j].fxMP[nth] -= force_bond_harmonic * dx;
		ions[j].fyMP[nth] -= force_bond_harmonic * dy;
		ions[j].fzMP[nth] -= force_bond_harmonic * dz;
		if(flags->eflag) vars->U_MP[nth].Uion+=rk*dr;
	}
}

void
PotentialAMBER::computeAngle(Variables *vars, FLAG *flags) {
	Atom *ions = vars->effectiveIn[0][0].inAtoms.data();
	Angle_type *ctypes = vars->ctypes.data();
	Angle *angles=vars->effectiveIn[0][0].angles.data();
	int asize=vars->effectiveIn[0][0].angles.size();
	#pragma omp parallel for
	for (int ian=0;ian<asize;ian++) {
		double dx1, dy1, dz1, dx2, dy2, dz2, rsq1, rsq2, r1, r2, C, Cs, dtheta, tk, a, a11, a12, a22, f1[3], f3[3];
		int nth=omp_get_thread_num();
		int i=angles[ian].atom1;
		int j=angles[ian].atom2;
		int k=angles[ian].atom3;
		int type=angles[ian].type;
		dx1 = ions[i].qx - ions[j].qx;
		dy1 = ions[i].qy - ions[j].qy;
		dz1 = ions[i].qz - ions[j].qz;
		rsq1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
		r1 = sqrt(rsq1);
		dx2 = ions[k].qx - ions[j].qx;
		dy2 = ions[k].qy - ions[j].qy;
		dz2 = ions[k].qz - ions[j].qz;
		rsq2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
		r2 = sqrt(rsq2);
		C = dx1*dx2 + dy1*dy2 + dz1*dz2;
		C /= r1*r2;
		Cs = 1/(sqrt(1.0-C*C));
		dtheta = acos(C) - ctypes[type].coeff[1];
		tk = ctypes[type].coeff[0] * dtheta;
		a = -2.0 * tk * Cs;
  	a11 = a*C / rsq1;
  	a12 = -a / (r1*r2);
  	a22 = a*C / rsq2;
		f1[0] = a11*dx1 + a12*dx2;
		f1[1] = a11*dy1 + a12*dy2;
  	f1[2] = a11*dz1 + a12*dz2;
  	f3[0] = a22*dx2 + a12*dx1;
  	f3[1] = a22*dy2 + a12*dy1;
  	f3[2] = a22*dz2 + a12*dz1;
  	ions[i].fxMP[nth] += f1[0];
		ions[i].fyMP[nth] += f1[1];
		ions[i].fzMP[nth] += f1[2];
  	ions[j].fxMP[nth] -= f1[0] + f3[0];
		ions[j].fyMP[nth] -= f1[1] + f3[1];
		ions[j].fzMP[nth] -= f1[2] + f3[2];
		ions[k].fxMP[nth] += f3[0];
		ions[k].fyMP[nth] += f3[1];
		ions[k].fzMP[nth] += f3[2];
    if (flags->eflag) vars->U_MP[nth].Uion+= tk*dtheta;
	}
}

void
PotentialAMBER::computeDihedral(Variables *vars, FLAG *flags) {
	Atom *ions = vars->effectiveIn[0][0].inAtoms.data();
	Dihedral_type *dtypes = vars->dtypes.data();
	Dihedral *dihedrals=vars->effectiveIn[0][0].dihedrals.data();
	int dsize=vars->effectiveIn[0][0].dihedrals.size();
	#pragma omp parallel for
	for (int idi=0;idi<dsize;idi++) {
		double ff2[3],ff4[3],ff1[3],ff3[3];

		int nth=omp_get_thread_num();
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
			if (flags->eflag) vars->U_MP[nth].Uion+= dtypes[type].coeff[JJ5] * p_;
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

		ions[i].fxMP[nth] += ff1[0];
		ions[i].fyMP[nth] += ff1[1];
		ions[i].fzMP[nth] += ff1[2];
		ions[j].fxMP[nth] += ff2[0];
		ions[j].fyMP[nth] += ff2[1];
		ions[j].fzMP[nth] += ff2[2];
		ions[k].fxMP[nth] += ff3[0];
		ions[k].fyMP[nth] += ff3[1];
		ions[k].fzMP[nth] += ff3[2];
		ions[l].fxMP[nth] += ff4[0];
		ions[l].fyMP[nth] += ff4[1];
		ions[l].fzMP[nth] += ff4[2];

	}
}

void
PotentialAMBER::initialAMBER(Variables *vars, FLAG *flags){
	std::vector<Pair> noLong;
	Pair p;
	for (auto &b : vars-> effectiveIn[0][0].bonds) {
		p.i=b.atom1;
		p.j=b.atom2;
		if(p.j<p.i) {
			p.i=b.atom2;
			p.j=b.atom1;
		}
		noLong.push_back(p);
	}
	for (auto &b : vars-> effectiveIn[0][0].angles) {
		p.i=b.atom1;
		p.j=b.atom3;
		if(p.j<p.i) {
			p.i=b.atom3;
			p.j=b.atom1;
		}
		noLong.push_back(p);
	}
	for (auto &b : vars-> effectiveIn[0][0].dihedrals) {
		p.i=b.atom1;
		p.j=b.atom4;
		if(p.j<p.i) {
			p.i=b.atom4;
			p.j=b.atom1;
		}
		noLong.push_back(p);
	}

	Atom *ions=vars->effectiveIn[0][0].inAtoms.data();
	const int is=vars->effectiveIn[0][0].inAtoms.size();
	for(int i=0;i<is-1;i++){
		for(int j=i+1;j<is;j++){
			int flag=1;
			for(auto &a : noLong){
				if(a.i==i && a.j==j){
					flag=0;
					break;
				}
			}
			if(flag==1) {
				p.i=i;
				p.j=j;
				longPair.push_back(p);
			}
		}
	}
}
