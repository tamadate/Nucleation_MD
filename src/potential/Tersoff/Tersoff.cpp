//------------------------------------------------------------------------
#include "potentialTersoff.hpp"
//------------------------------------------------------------------------



/////////////////////////////////////////////////////////////////////
/*
	- Calculate the intra-atomic interaction (Stillinger-Weber)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialTersoff::compute(Variables *vars, FLAG *flags) {
	Atom *ions = vars->Molecules[0].inAtoms.data();
	const int is = vars->Molecules[0].inAtoms.size();
	U=0;
	Pressure=0;
	for(auto &a : pairs){
		int i=a.i;
		neighshort.clear();
//	repulsive forces d(fCfR)/dr
		for(auto j : a.j){
			delr[0]=ions[i].qx-ions[j].qx;
			delr[1]=ions[i].qy-ions[j].qy;
			delr[2]=ions[i].qz-ions[j].qz;
		//	adjust_periodic(delr[0], delr[1], delr[2]);
			double rsq = (delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2]);
			if(rsq>cut2) continue;
			double r=sqrt(rsq);
			double force_pair=-(fC(r)*fR_d(r)+fR(r)*fC_d(r))/r*0.5;
			ions[i].fx+=force_pair*delr[0];
			ions[i].fy+=force_pair*delr[1];
			ions[i].fz+=force_pair*delr[2];
			ions[j].fx-=force_pair*delr[0];
			ions[j].fy-=force_pair*delr[1];
			ions[j].fz-=force_pair*delr[2];

			neighshort.push_back(j);
			U+=fC(r)*fR(r)*0.5;
		}

//	three body interaction
		for(auto j : neighshort){
			delr[0]=ions[i].qx-ions[j].qx;
			delr[1]=ions[i].qy-ions[j].qy;
			delr[2]=ions[i].qz-ions[j].qz;
		//	adjust_periodic(delr[0], delr[1], delr[2]);
			double rsq = (delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2]);
			double r=sqrt(rsq);
			double zeta_ij=0;

			for (auto k : neighshort){
				if(k==j) continue;
				double dx2=ions[i].qx-ions[k].qx;
				double dy2=ions[i].qy-ions[k].qy;
				double dz2=ions[i].qz-ions[k].qz;
			//	adjust_periodic(dx2, dy2, dz2);
				double rsq2=(dx2*dx2+dy2*dy2+dz2*dz2);
				if(rsq2>cut2) continue;
				double r2=sqrt(rsq2);
				zeta_ij+=zeta(r,r2,delr[0]*dx2,delr[1]*dy2,delr[2]*dz2);
			}
			double prefactor=0;
			double force_pair=-0.5*force_zeta(r,zeta_ij,prefactor);
			ions[i].fx += force_pair * delr[0];
			ions[i].fy += force_pair * delr[1];
			ions[i].fz += force_pair * delr[2];
			ions[j].fx -= force_pair * delr[0];
			ions[j].fy -= force_pair * delr[1];
			ions[j].fz -= force_pair * delr[2];

			U+=bij(zeta_ij)*fC(r)*fA(r)*0.5;

			if(zeta_ij==0) continue;
			for (auto k : neighshort) {
				if(k==j) continue;
				delr2[0] = ions[i].qx - ions[k].qx;
				delr2[1] = ions[i].qy - ions[k].qy;
				delr2[2] = ions[i].qz - ions[k].qz;
			//	adjust_periodic(delr2[0], delr2[1], delr2[2]);
				double rsq2 = (delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2]);
				if(rsq2>cut2) continue;
				double r2 = sqrt(rsq2);
				double fi[3],fj[3],fk[3];
				attractive(prefactor,r,r2,delr,delr2,fi,fj,fk);
				ions[i].fx+=fi[0]*0.5;
				ions[i].fy+=fi[1]*0.5;
				ions[i].fz+=fi[2]*0.5;
				ions[j].fx+=fj[0]*0.5;
				ions[j].fy+=fj[1]*0.5;
				ions[j].fz+=fj[2]*0.5;
				ions[k].fx+=fk[0]*0.5;
				ions[k].fy+=fk[1]*0.5;
				ions[k].fz+=fk[2]*0.5;
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
}

double
PotentialTersoff::compute_virial(Variables *vars) {
/*	Atom *ions = vars->Molecules[0].inAtoms.data();
	const int is = vars->Molecules[0].inAtoms.size();
	double virial=0;
	for(auto &a : pairs){
		int i=a.i;
		int j=a.j;
	repulsive forces d(fCfR)/dr
		double delr[3];
		delr[0]=ions[i].qx-ions[j].qx;
		delr[1]=ions[i].qy-ions[j].qy;
		delr[2]=ions[i].qz-ions[j].qz;
		adjust_periodic(delr[0], delr[1], delr[2]);
		double rsq = (delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2]);
		if(rsq>cut2) continue;
		double r=sqrt(rsq);
		double force_pair=-(fC(r)*fR_d(r)+fR(r)*fC_d(r));
		virial+=force_pair*r;

	three body interaction
		double zeta_ij=0;
		for (auto &k : a.k){
			double dx2=ions[i].qx-ions[k].qx;
			double dy2=ions[i].qy-ions[k].qy;
			double dz2=ions[i].qz-ions[k].qz;
			adjust_periodic(dx2, dy2, dz2);
			double rsq2=(dx2*dx2+dy2*dy2+dz2*dz2);
			if(rsq2>cut2) continue;
			double r2=sqrt(rsq2);
			zeta_ij+=zeta(r,r2,delr[0]*dx2,delr[1]*dy2,delr[2]*dz2);
		}
		double prefactor = -1*fA(r)*fC(r)*bij_d(zeta_ij);
		if (zeta_ij==0) prefactor=0;
		force_pair=-1*bij(zeta_ij)*(fA(r)*fC_d(r)+fA_d(r)*fC(r));
		virial+=force_pair*r;

		if(zeta_ij==0) continue;
		for (auto &k : a.k){
			double delr2[3];
			delr2[0] = ions[i].qx - ions[k].qx;
			delr2[1] = ions[i].qy - ions[k].qy;
			delr2[2] = ions[i].qz - ions[k].qz;
			adjust_periodic(delr2[0], delr2[1], delr2[2]);
			double rsq2 = (delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2]);
			if(rsq2>cut2) continue;
			double r2 = sqrt(rsq2);
			double rij_hat[3],rik_hat[3];
			for (int i=0; i<3; i++) {
				rij_hat[i]=delr[i]/r;
				rik_hat[i]=delr2[i]/r2;
			}
			double fc=fC(r2);
			double dfc=fC_d(r2);
			double cos_theta=rij_hat[0]*rik_hat[0]+rij_hat[1]*rik_hat[1]+rij_hat[2]*rik_hat[2];
			double g = gijk(cos_theta);
			double g_d = gijk_d(cos_theta);
			virial-=g*dfc*r2*prefactor;
			virial-=fc*g_d*(r2-cos_theta*r)/r*prefactor;
			virial-=fc*g_d*(r-cos_theta*r2)/r2*prefactor;
		}
	}
	return virial*23.06;*/
	return 0;
}



void
PotentialTersoff::attractive(double prefactor, double rij, double rik, double *delrij, double *delrik, double *fi, double *fj, double *fk){
	double rij_hat[3],rik_hat[3];
	for (int i=0; i<3; i++) {
		rij_hat[i]=delrij[i]/rij;
		rik_hat[i]=delrik[i]/rik;
	}
	zeta_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk);
}

void
PotentialTersoff::zeta_d(double prefactor, double *rij_hat, double rij,double *rik_hat, double rik,double *dri, double *drj, double *drk){
	double g,g_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
	double dcosdri[3],dcosdrj[3],dcosdrk[3];
	fc=fC(rik);
	dfc=fC_d(rik);
	cos_theta=rij_hat[0]*rik_hat[0]+rij_hat[1]*rik_hat[1]+rij_hat[2]*rik_hat[2];
	g = gijk(cos_theta);
	g_d = gijk_d(cos_theta);
	for (int i=0;i<3;i++) {
		dcosdrj[i]=(rik_hat[i]-cos_theta*rij_hat[i])/rij;
		dcosdrk[i]=(rij_hat[i]-cos_theta*rik_hat[i])/rik;
		dcosdri[i]=-dcosdrj[i]-dcosdrk[i];
	}
	for (int i=0;i<3;i++) {
		dri[i]=-dfc*g*rik_hat[i];
		dri[i]+=fc*g_d*dcosdri[i];
		dri[i]*=prefactor;
	}
	for (int i=0;i<3;i++) {
		drj[i]=fc*g_d*dcosdrj[i];
		drj[i]*=prefactor;
	}
	for (int i=0;i<3;i++) {
		drk[i]=dfc*g*rik_hat[i];
		drk[i]+=fc*g_d*dcosdrk[i];
		drk[i]*=prefactor;
	}
}

double
PotentialTersoff::force_zeta(double r, double zeta_ij, double &prefactor){
	prefactor = fA(r)*fC(r)*bij_d(zeta_ij);
	if (zeta_ij==0) prefactor=0;
	return bij(zeta_ij)*(fA(r)*fC_d(r)+fA_d(r)*fC(r))/r;
}

double
PotentialTersoff::zeta(double rij, double rik, double x1x2, double y1y2, double z1z2) {
	double costheta=(x1x2+y1y2+z1z2)/(rij*rik);
	double arg=exp(lamda3*lamda3*lamda3*pow(rij-rik,3));
	return arg*gijk(costheta)*fC(rik);
}

double
PotentialTersoff::gijk(double costheta) {
	double hcth= costheta0-costheta;
	return gamma*(1.0+c2/d2-c2/(d2+hcth*hcth));
}

double
PotentialTersoff::gijk_d(double costheta) {
	double hcth= costheta0-costheta;
	double numerator=-2.0*c2*hcth;
	double denominator=1/(d2+hcth*hcth);
	return gamma*numerator*denominator*denominator;
}

double
PotentialTersoff::bij(double zeta) {

	double tmp = beta * zeta;
/*	if (tmp > ter_c1) return 1.0/sqrt(tmp);
	if (tmp > ter_c2) return (1.0 - pow(tmp,-ter_n) / (2.0*ter_n))/sqrt(tmp);
	if (tmp < ter_c4) return 1.0;
	if (tmp < ter_c3) return 1.0 - pow(tmp,ter_n)/(2.0*ter_n);*/
	if(zeta==0) return 1.0;
	else return pow(1.0 + pow(tmp,ter_n), -1.0/(2.0*ter_n));
}

double
PotentialTersoff::bij_d(double zeta){
	double tmp = beta * zeta;
/*	if (tmp > ter_c1) return beta*-0.5*pow(tmp,-1.5);
	if (tmp > ter_c2) return beta*(-0.5*pow(tmp,-1.5)*(1.0-(1.0+1.0/(2.0*ter_n))*pow(tmp,-ter_n)));
	if (tmp < ter_c4) return 0.0;
	if (tmp < ter_c3) return -0.5*beta * pow(tmp,ter_n-1.0);*/
	double tmp_n = pow(tmp,ter_n);
	return -0.5*pow(1.0+tmp_n, -1.0-(1.0/(2.0*ter_n)))*tmp_n / zeta;
}

void
PotentialTersoff::check_pairlist(Variables *vars){
	loop_t++;
	if(loop_t>loop_update_t) make_pair(vars);
}

void
PotentialTersoff::make_pair(Variables *vars){
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
		//	adjust_periodic(dx, dy, dz);
			double rsq = (dx * dx + dy * dy + dz * dz);
			if (rsq < T_ML2) js.push_back(j);
		}
		p.j=js;
		pairs.push_back(p);
	}
	loop_t=0;
}
