//------------------------------------------------------------------------
#include "potentialSW.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*
	- Calculate the intra-atomic interaction (Stillinger-Weber)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialSW::compute(Variables *vars, FLAG *flags) {
	Atom *ions = vars->Molecules[0].inAtoms.data();
	const int is = vars->Molecules[0].inAtoms.size();
	U=0;
	Pressure=0;
	virial=0;

	for(auto &a : pairs){
		int i=a.i;
		neighshort.clear();
//	repulsive forces d(fCfR)/dr
		for(auto j : a.j){
			delr[0]=ions[i].qx-ions[j].qx;
			delr[1]=ions[i].qy-ions[j].qy;
			delr[2]=ions[i].qz-ions[j].qz;
			//adjust_periodic(delr[0], delr[1], delr[2]);
			double rsq = (delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2]);
			if(rsq>cut2) continue;
			double force_pair=twobody(rsq)*0.5;
			ions[i].fx+=force_pair*delr[0];
			ions[i].fy+=force_pair*delr[1];
			ions[i].fz+=force_pair*delr[2];
			ions[j].fx-=force_pair*delr[0];
			ions[j].fy-=force_pair*delr[1];
			ions[j].fz-=force_pair*delr[2];

			neighshort.push_back(j);

//			virial+=force_pair*r*r;
		}

//	three body interaction
		for(auto j : neighshort){
			delr[0]=ions[i].qx-ions[j].qx;
			delr[1]=ions[i].qy-ions[j].qy;
			delr[2]=ions[i].qz-ions[j].qz;
			//adjust_periodic(delr[0], delr[1], delr[2]);
			double rsq = (delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2]);

//			virial+=force_pair*r*r;

			for (auto k : neighshort){
				if(k==j) continue;
				delr2[0] = ions[i].qx - ions[k].qx;
				delr2[1] = ions[i].qy - ions[k].qy;
				delr2[2] = ions[i].qz - ions[k].qz;
			//	adjust_periodic(delr2[0], delr2[1], delr2[2]);
				double rsq2 = (delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2]);
				if(rsq2>cut2) continue;
				double fi[3],fj[3],fk[3];
				threebody(rsq,rsq2,delr,delr2,fj,fk);
				ions[i].fx+=(fj[0]+fk[0])*0.5;
				ions[i].fy+=(fj[1]+fk[1])*0.5;
				ions[i].fz+=(fj[2]+fk[2])*0.5;
				ions[j].fx-=fj[0]*0.5;
				ions[j].fy-=fj[1]*0.5;
				ions[j].fz-=fj[2]*0.5;
				ions[k].fx-=fk[0]*0.5;
				ions[k].fy-=fk[1]*0.5;
				ions[k].fz-=fk[2]*0.5;
			}
		}
	}

	for(int i=0;i<is;i++){
/*	FILE*f=fopen("force.dat", "a");
	fprintf(f, "%f\t%f\t%f\t%e\t%e\t%e\t%f\n", ions[i].qx, ions[i].qy, ions[i].qz, ions[i].fx, ions[i].fy, ions[i].fz, energy);
	fclose(f);*/
		ions[i].fx*=23.06;
		ions[i].fy*=23.06;
		ions[i].fz*=23.06;
	}
	U*=0.5;
//	vars->Potential+=energy*23.06;
//	vars->virial+=virial*23.06;

}

double
PotentialSW::computeVirial(Variables *vars) {
	Atom *ions = vars->Molecules[0].inAtoms.data();
	const int is = vars->Molecules[0].inAtoms.size();
/*	double virial=0;
	for(auto &a : vars->ions){
		virial+=a.fx*a.qx;
		virial+=a.fy*a.qy;
		virial+=a.fz*a.qz;
	}*/
	return virial*23.061;
}



double
PotentialSW::twobody(double rsq){
	double r,rinvsq,rp,rq,rainv,rainvsq,expsrainv;

	r = sqrt(rsq);
	rinvsq = 1.0/rsq;
	rp = pow(r,-powerp);
	rq = pow(r,-powerq);
	rainv = 1.0 / (r - cut);
	rainvsq = rainv*rainv*r;
	expsrainv = exp(sigma * rainv);
	U+=(c5*rp-c6*rq)*expsrainv;
	return (c1*rp - c2*rq + (c3*rp -c4*rq) * rainvsq) * expsrainv * rinvsq;
}


void
PotentialSW::threebody(double rsq1, double rsq2,double *rij, double *rik,double *fj, double *fk)
{
	double r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
	double r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
	double rinv12,cs,delcs,delcssq,facexp,facrad,frad1,frad2;
	double facang,facang12,csfacang,csfac1,csfac2;

	r1 = sqrt(rsq1);
	rinvsq1 = 1.0/rsq1;
	rainv1 = 1.0/(r1 - cut);
	gsrainv1 = sigma_gamma * rainv1;
	gsrainvsq1 = gsrainv1*rainv1/r1;
	expgsrainv1 = exp(gsrainv1);

	r2 = sqrt(rsq2);
	rinvsq2 = 1.0/rsq2;
	rainv2 = 1.0/(r2 - cut);
	gsrainv2 = sigma_gamma * rainv2;
	gsrainvsq2 = gsrainv2*rainv2/r2;
	expgsrainv2 = exp(gsrainv2);

	rinv12 = 1.0/(r1*r2);
	cs = (rij[0]*rik[0] + rij[1]*rik[1] + rij[2]*rik[2]) * rinv12;
	delcs = cs - costheta0;
	delcssq = delcs*delcs;

	facexp = expgsrainv1*expgsrainv2;

	facrad = lambda_epsilon * facexp*delcssq;
	frad1 = facrad*gsrainvsq1;
	frad2 = facrad*gsrainvsq2;
	facang = lambda_epsilon2 * facexp*delcs;
	facang12 = rinv12*facang;
	csfacang = cs*facang;

	csfac1 = rinvsq1*csfacang;

	fj[0] = rij[0]*(frad1+csfac1)-rik[0]*facang12;
	fj[1] = rij[1]*(frad1+csfac1)-rik[1]*facang12;
	fj[2] = rij[2]*(frad1+csfac1)-rik[2]*facang12;

	csfac2 = rinvsq2*csfacang;

	fk[0] = rik[0]*(frad2+csfac2)-rij[0]*facang12;
	fk[1] = rik[1]*(frad2+csfac2)-rij[1]*facang12;
	fk[2] = rik[2]*(frad2+csfac2)-rij[2]*facang12;

	U+=facrad;
}

void
PotentialSW::check_pairlist(Variables *vars){
	loop_t++;
	if(loop_t>loop_update_t) make_pair(vars);
}

void
PotentialSW::make_pair(Variables *vars){
	pairs.clear();
	Atom *ions = vars->Molecules[0].inAtoms.data();
	int is=vars->Molecules[0].inAtoms.size();
	for (int i=0; i<is; i++){
		Pair_many p;
		p.i=i;
		vector<int> js;
		for (int j=0; j<is; j++){
			if (i==j) continue;
			double dx = ions[i].qx - ions[j].qx;
			double dy = ions[i].qy - ions[j].qy;
			double dz = ions[i].qz - ions[j].qz;
	//		adjust_periodic(dx, dy, dz);
			double rsq = (dx * dx + dy * dy + dz * dz);
			if (rsq < sw_ML2) js.push_back(j);
		}
		p.j=js;
		pairs.push_back(p);
	}
	loop_t=0;
}
