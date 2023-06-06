#include "../potential.hpp"


/*########################################################################################

-----compute intramolecular interaction-----

#######################################################################################*/

/**********************************Force calculation******************************************/
void
PotentialBorn::compute(Variables *vars) {
	Atom *ions = vars->ions.data();
	const int is = vars->ions.size();
	vars->times.tion-=omp_get_wtime();
	for(int i=0; i<is-1; i++){
		for(int j=i+1; j<is; j++){
			double dx = ions[i].qx - ions[j].qx;
			double dy = ions[i].qy - ions[j].qy;
			double dz = ions[i].qz - ions[j].qz;
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			double r6inv = r2inv * r2inv * r2inv;
			double r=sqrt(rsq);
			int type1=ions[i].type;
			int type2=ions[j].type;

//vars->bornCoeff[type1][type2][0]=A
//vars->bornCoeff[type1][type2][1]=6C
//vars->bornCoeff[type1][type2][2]=8D
//vars->bornCoeff[type1][type2][3]=sigma
//vars->bornCoeff[type1][type2][4]=1/rho

			double rexp=exp(bornCoeff[type1][type2][4]*(bornCoeff[type1][type2][3]-r));
			double force1 = bornCoeff[type1][type2][4]*bornCoeff[type1][type2][0]*r*rexp;
			double force2 = -bornCoeff[type1][type2][1]*r6inv;
			double force3 = -bornCoeff[type1][type2][2]*r6inv*r2inv;
			double force_coul = qqrd2e*ions[i].charge*ions[j].charge*sqrt(r2inv);
			double force_pair = (force1+force2+force3+force_coul)*r2inv;

			ions[i].fx += force_pair * dx;
			ions[i].fy += force_pair * dy;
			ions[i].fz += force_pair * dz;
			ions[j].fx -= force_pair * dx;
			ions[j].fy -= force_pair * dy;
			ions[j].fz -= force_pair * dz;
			if(vars->eflag) {
				vars->U.Uion+=force_coul;
				vars->U.Uion+=rexp*bornCoeff[type1][type2][0];
				vars->U.Uion-=bornCoeff[type1][type2][1]/6.0*r6inv;
				vars->U.Uion-=bornCoeff[type1][type2][2]/8.0*r6inv*r2inv;
			}
			//vars->totalVirial+=force_lj;
		}
	}
	vars->times.tion+=omp_get_wtime();
}


PotentialBorn::PotentialBorn(Variables *vars){
	//difine Born-Mayer-Huggins conefficients for "only" NaCl
	//https://doi.org/10.1063/1.1522375
	//https://doi.org/10.1016/0022-3697(64)90160-X
	//0:A/rho, 1:6C, 2:8D, 3:sigma, 4:1/rho
	//NaNa, A=25.4435kJ/mol, C=101.1719kJ/mol, D=48.1771kJ/mol, sigma=2.340A, 1/rho=3.1546A-1
	//NaCl, A=20.3548kJ/mol, C=674.4793kJ/mol, D=837.077kJ/mol, sigma=2.755A, 1/rho=3.1546A-1
	//ClCl, A=15.2661kJ/mol, C=6985.6786kJ/mol, D=14031.5785kJ/mol, sigma=3.170A, 1/rho=3.1546A-1

  	bornCoeff[0][0][0]=25.4435*3.1546/4.184;
  	bornCoeff[0][1][0]=bornCoeff[1][0][0]=20.3548*3.1546/4.184;
  	bornCoeff[1][1][0]=15.2661*3.1546/4.184;

  	bornCoeff[0][0][1]=6*101.1719/4.184;
  	bornCoeff[0][1][1]=bornCoeff[1][0][1]=6*674.4793/4.184;
  	bornCoeff[1][1][1]=6*6985.6786/4.184;

  	bornCoeff[0][0][2]=8*48.1771/4.184;
  	bornCoeff[0][1][2]=bornCoeff[1][0][2]=8*837.077/4.184;
  	bornCoeff[1][1][2]=8*14031.5785/4.184;

  	bornCoeff[0][0][3]=2.340;
  	bornCoeff[0][1][3]=bornCoeff[1][0][3]=2.755;
  	bornCoeff[1][1][3]=3.170;

  	bornCoeff[0][0][4]=3.1546;
  	bornCoeff[0][1][4]=bornCoeff[1][0][4]=3.1546;
  	bornCoeff[1][1][4]=3.1546;
}