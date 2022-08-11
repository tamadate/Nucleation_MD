#pragma once
#include "variables.hpp"
#include "flags.hpp"

//------------------------------------------------------------------------
class Potential {
	private:
		string potName="potential";
	public:
		virtual void printName(void) {cout<<potName<<endl;}
		virtual void compute(Variables *vars, FLAG *flags);
		Potential(){};
		~Potential(){};
};

class PotentialGasIon : public Potential {
	private:
		string potName="LJ gas-ion";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialGasIon(){};
		~PotentialGasIon(){};
};

class PotentialVaporVapor : public Potential {
	private:
		string potName="LJ vapor-vapor";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialVaporVapor(){};
		~PotentialVaporVapor(){};	
};

class PotentialVaporGas : public Potential {
	private:
		string potName="LJ vapor-gas";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialVaporGas(){};
		~PotentialVaporGas(){};	
};

class PotentialVaporIon : public Potential {
	private:
		string potName="LJ vapor-ion";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialVaporIon(){};
		~PotentialVaporIon(){};	
};

class PotentialGasGas : public Potential {
	private:
		string potName="LJ gas-gas";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialGasGas(){};
		~PotentialGasGas(){};	
};


class PotentialAMBER : public Potential {
	private:
		string potName="AMBER ion";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		void computeBond(Variables *vars, FLAG *flags);
		void computeAngle(Variables *vars, FLAG *flags);
		void computeDihedral(Variables *vars, FLAG *flags);
		PotentialAMBER(){};
		~PotentialAMBER(){};	
};


class PotentialGasIntra : public Potential {
	private:
		string potName="Gas intra";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialGasIntra(){};
		~PotentialGasIntra(){};	
};

class PotentialBorn : public Potential {
	private:
		string potName="Ion BMH";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialBorn(){};
		~PotentialBorn(){};	
};

class PotentialEfield : public Potential {
	private:
		string potName="Ion E-filed";
	public:
		void printName(void) {cout<<potName<<endl;}
		double Ecoeff[3];
		void compute(Variables *vars, FLAG *flags);
		PotentialEfield(){};
		~PotentialEfield(){};	
};

class PotentialIonDipole : public Potential {
	private:
		string potName="Induced dipole ion-gas";
	public:
		void printName(void) {cout<<potName<<endl;}
		double alphagas;
		double zion;
		void compute(Variables *vars, FLAG *flags);
		PotentialIonDipole(){};
		~PotentialIonDipole(){};	
};

class PotentialVaporIntra : public Potential {
	private:
		string potName="AMBER vapor";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialVaporIntra(){};
		~PotentialVaporIntra(){};	
};

class PotentialIntraTIP3P : public Potential {
	private:
		string potName="TIP3P vapor";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialIntraTIP3P(){};
		~PotentialIntraTIP3P(){};	
};



class PotentialSW : public Potential {
	private:
		double computeVirial(Variables *vars);
		double twobody(double rsq);
		void threebody(double rsq1, double rsq2,double *r1, double *r2,double *fj, double *fk);
		void check_pairlist(Variables *vars);
		void make_pair(Variables *vars);

		const double epsilon=2.1683;
		const double sigma=2.0951;
		const double a=1.80;
		const double lambda=21.0;
		const double gamma=1.20;
		const double costheta0=-0.333333333333;
		const double A=7.049556277;
		const double B=0.6022245584;
		const double powerp=4.0;
		const double powerq=0.0;
		const double tol=0.0;
		const double c1=A*epsilon*powerp*B*pow(sigma,powerp);
		const double c2=A*epsilon*powerq*pow(sigma,powerq);
		const double c3=A*epsilon*B*pow(sigma,powerp+1);
		const double c4=A*epsilon*pow(sigma,powerq+1);
		const double c5=A*epsilon*B*pow(sigma,powerp);
		const double c6=A*epsilon*pow(sigma,powerq);
		const double sigma_gamma=sigma*gamma;
		const double lambda_epsilon=lambda*epsilon;
		const double lambda_epsilon2=2*lambda*epsilon;

		const double cut=a*sigma;
		const double cut2=cut*cut;
		const double MARGIN=5;
		const double sw_ML2=(cut+MARGIN)*(cut+MARGIN);

		std::vector<Pair> pairs_inner;

		double delr[3],delr2[3];

		int loop_t;
		const int loop_update_t=100;
		std::vector<Pair_many> pairs;
		std::vector<int> neighshort;

	public:
		void compute(Variables *vars, FLAG *flags);
		double energy;
		double Pressure;
		double virial;


		PotentialSW(){};
		~PotentialSW(){};	
};


class PotentialTersoff : public Potential{
	public:
		void compute(Variables *vars, FLAG *flags);
		double compute_virial(Variables *vars);
		void check_pairlist(Variables *vars);
		void make_pair(Variables *vars);
		double energy;
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
