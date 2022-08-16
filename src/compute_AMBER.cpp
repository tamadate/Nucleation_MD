#include "potential.hpp"
/*########################################################################################

-----compute intramolecular interaction-----

#######################################################################################*/

/**********************************Force calculation******************************************/
void
PotentialAMBER::compute(Variables *vars, FLAG *flags) {
	computeLong(vars,flags);
	cout<<"bond"<<endl;
	computeBond(vars,flags);
	cout<<"angle"<<endl;
	computeAngle(vars,flags);
	cout<<"dihedral"<<endl;
	computeDihedral(vars,flags);
	cout<<"done"<<endl;
}


void
PotentialAMBER::computeLong(Variables *vars, FLAG *flags) {
	Atom *ions = vars->ions.data();
/*intra-molecular interaction (bond)*/
	Bond_type *btypes = vars->btypes.data();
	for (auto &p : longPair) {
		int i=p.i, j=p.j;
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
		ions[i].fx += force_pair * dx;
		ions[i].fy += force_pair * dy;
		ions[i].fz += force_pair * dz;
		ions[j].fx -= force_pair * dx;
		ions[j].fy -= force_pair * dy;
		ions[j].fz -= force_pair * dz;
		if(flags->eflag) {
			vars->Uion+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			vars->Uion+=force_coul;
		}
	}
}

void
PotentialAMBER::computeAngle(Variables *vars, FLAG *flags) {
	Atom *ions = vars->ions.data();
/*intra-molecular interaction (angle)*/
	double dx1, dy1, dz1, dx2, dy2, dz2, rsq1, rsq2, r1, r2, C, Cs, dtheta, tk, a, a11, a12, a22, f1[3], f3[3];
	Angle_type *ctypes = vars->ctypes.data();
	for (auto &c : vars-> angles) {
        int i, j, k, type;
		i=c.atom1, j=c.atom2, k=c.atom3, type=c.type;
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
    	ions[i].fx += f1[0];
		ions[i].fy += f1[1];
		ions[i].fz += f1[2];
    	ions[j].fx -= f1[0] + f3[0];
		ions[j].fy -= f1[1] + f3[1];
		ions[j].fz -= f1[2] + f3[2];
		ions[k].fx += f3[0];
		ions[k].fy += f3[1];
		ions[k].fz += f3[2];
	    if (flags->eflag) vars->Uion+= tk*dtheta;
	}
}

void
PotentialAMBER::computeDihedral(Variables *vars, FLAG *flags) {
	Atom *ions = vars->ions.data();
/*intra-molecular interaction (dihedral)*/
	double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
	double edihedral,ff2[3],ff4[3],ff1[3],ff3[3];
	double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
	double df,df1,ddf1,fg,hg,fga,hgb,gaa,gbb;
	double dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;
	double c,s,p_,sx2,sy2,sz2, m;
	Dihedral_type *dtypes = vars->dtypes.data();
	int loop=0;
	for (auto &d : vars-> dihedrals) {
		loop++;
		cout<<"dihedral"<<loop<<endl;
		int i=d.atom1, j=d.atom2, k=d.atom3, l=d.atom4, type=d.type;
   //     cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<endl;
		// 1st bond
		vb1x = ions[i].qx - ions[j].qx;
		vb1y = ions[i].qy - ions[j].qy;
		vb1z = ions[i].qz - ions[j].qz;
		// 2nd bond
		vb2x = ions[k].qx - ions[j].qx;
		vb2y = ions[k].qy - ions[j].qy;
		vb2z = ions[k].qz - ions[j].qz;
		vb2xm = -vb2x;
		vb2ym = -vb2y;
		vb2zm = -vb2z;
		// 3rd bond
		vb3x = ions[l].qx - ions[k].qx;
		vb3y = ions[l].qy - ions[k].qy;
		vb3z = ions[l].qz - ions[k].qz;

		ax = vb1y*vb2zm - vb1z*vb2ym;
		ay = vb1z*vb2xm - vb1x*vb2zm;
		az = vb1x*vb2ym - vb1y*vb2xm;
		bx = vb3y*vb2zm - vb3z*vb2ym;
		by = vb3z*vb2xm - vb3x*vb2zm;
		bz = vb3x*vb2ym - vb3y*vb2xm;
		rasq = ax*ax + ay*ay + az*az;
		rbsq = bx*bx + by*by + bz*bz;
		rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
		rg = sqrt(rgsq);
		rginv = ra2inv = rb2inv = 0.0;
		if (rg > 0) rginv = 1.0/rg;
		if (rasq > 0) ra2inv = 1.0/rasq;
		if (rbsq > 0) rb2inv = 1.0/rbsq;
		rabinv = sqrt(ra2inv*rb2inv);
		c = (ax*bx + ay*by + az*bz)*rabinv;
		s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

		df = 0.0;
		int J= dtypes[type].multi;
		cout<<J<<endl;

		for(int JJ=0;JJ<J;JJ++){
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
			if (flags->eflag) vars->Uion+= dtypes[type].coeff[JJ5] * p_;
		}

       // cout<<df<<endl;
		fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
		hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
		fga = fg*ra2inv*rginv;
		hgb = hg*rb2inv*rginv;
		gaa = -ra2inv*rg;
		gbb = rb2inv*rg;
		dtfx = gaa*ax;
		dtfy = gaa*ay;
		dtfz = gaa*az;
		dtgx = fga*ax - hgb*bx;
		dtgy = fga*ay - hgb*by;
		dtgz = fga*az - hgb*bz;
		dthx = gbb*bx;
		dthy = gbb*by;
		dthz = gbb*bz;
		sx2 = df*dtgx;
		sy2 = df*dtgy;
		sz2 = df*dtgz;
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

void initialAMBER(Variables *vars, FLAG *flags){
	std::vector<Pair> noLong;
	Pair p;
	for (auto &b : vars-> bonds) {
		p.i=b.atom1;
		p.j=b.atom2;
		if(p.j<p.i) {
			p.i=b.atom2;
			p.j=b.atom1;
		}
		noLong.push_back(p)
	}
	for (auto &b : vars-> angles) {
		p.i=b.atom1;
		p.j=b.atom3;
		if(p.j<p.i) {
			p.i=b.atom3;
			p.j=b.atom1;
		}
		noLong.push_back(p)
	}
	for (auto &b : vars-> dihedrals) {
		p.i=b.atom1;
		p.j=b.atom4;
		if(p.j<p.i) {
			p.i=b.atom4;
			p.j=b.atom1;
		}
		noLong.push_back(p)
	}

	Atom *ions=vars->ions.data();
	const int is=vars->ions.size();
	for(int i=0;i<is-1;i++){
		for(int j=i+1;j<is;j++){
			p.i=i;
			p.i=j;
			int flag=1;
			for(auto &a : noLong){
				if(a.i==p.i && a.j==p.j){
					flag=0;
					break;
				}
			}
			if(flag) longPair.push_back(p);
		}
	}
}
