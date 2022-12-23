#include "../../variables.hpp"
#include "../basePotential.hpp"

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
		double U;
		double Pressure;
		double virial;


		PotentialSW(){};
		~PotentialSW(){};
};
