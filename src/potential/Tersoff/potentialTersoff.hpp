#include "../../variables.hpp"
#include "../basePotential.hpp"

class PotentialTersoff : public Potential{
	public:
		void compute(Variables *vars, FLAG *flags);
		double compute_virial(Variables *vars);
		void check_pairlist(Variables *vars);
		void make_pair(Variables *vars);
		double U;
		double Pressure;

	private:
		const double gamma=1.0;
		const double lamda3=0;
		const double c=100390;
		const double c2=c*c;
		const double d=16.218;
		const double d2=d*d;
		const double costheta0=-0.59826;
		const double ter_n=0.78734;
		const double ter_c1=pow(2.0*ter_n*1.0e-16,-1.0/ter_n);
		const double ter_c2=pow(2.0*ter_n*1.0e-8,-1.0/ter_n);
		const double ter_c3=1/ter_c2;
		const double ter_c4=1/ter_c1;
		const double beta=1.0999e-6;
		const double lamda2=1.7322;
		const double B=471.18;
		const double R=2.85;
		const double D=0.15;
		const double cut=R+D;
		const double cut2=cut*cut;
		const double cutl=R-D;
		const double cutl2=cutl*cutl;
		const double lamda1=2.4799;
		const double A=1830.8;
		const double MARGIN=5;
		const double T_ML2=(R+D+MARGIN)*(R+D+MARGIN);

		std::vector<Pair> pairs_inner;

		/*******************FUNCTION************************/
		double gijk_d(double costheta);
		double gijk(double costheta);
		double zeta(double rij, double rik, double x1x2, double y1y2, double z1z2);
	  	void zeta_d(double prefactor, double *rij_hat, double rij, double *rik_hat, double rik, double *fi, double *fj, double *fk);
		double bij(double zeta);
		double bij_d(double zeta);
		double fA(double r)	{return -B*exp(-lamda2*r);}
		double fA_d(double r)	{return B*lamda2*exp(-lamda2*r);}
		double fR(double r) {return A*exp(-lamda1*r);}
		double fR_d(double r) {return -A*lamda1*exp(-lamda1*r);}
		double fC(double r) {
			double fc;
			fc=0.5*(1-sin(M_PI*0.5*(r-R)/D));
			if (r<cutl) fc=1;
			if (r>cut) fc=0;
			return fc;
		}
		double fC_d(double r) {
			double fc;
			fc=-0.25*M_PI/D*cos(M_PI*0.5*(r-R)/D);
			if (r<cutl || r>cut) fc=0;
			return fc;
		}
		void attractive(double prefactor, double rij, double rik, double *delrij, double *delrik, double *fi, double *fj, double *fk);
		double force_zeta(double r, double zeta_ij, double &prefactor);
		double delr[3],delr2[3];

		int loop_t;
		const int loop_update_t=100;
		std::vector<Pair_many> pairs;
		std::vector<int> neighshort;

};


//------------------------------------------------------------------------
