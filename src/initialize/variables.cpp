#include "../variables.hpp"

Variables::Variables(void) {
  time = 0.0;
  U.Uion=U.Ugas=U.Uvap=U.Ugi=U.Ugg=U.Uvg=U.Uvi=U.Uvv=0;
  times.tion=times.tgas=times.tvap=times.tgi=times.tvv=times.tvg=times.tvi=times.tpair=0;
}


void
Variables::setCrossPotentials(void){
  int num_atoms=atypes.size();
	pair_coeff.resize(num_atoms);
	for (int i=0;i<num_atoms;i++){
		pair_coeff[i].resize(num_atoms);
		for (int j=0;j<num_atoms;j++){
			pair_coeff[i][j].resize(2);
		}
	}
	for (int i=0;i<num_atoms;i++){
		for (int j=0;j<num_atoms;j++){
			double epu=sqrt(atypes[i].coeff1*atypes[j].coeff1);
			double sigma=(atypes[i].coeff2+atypes[j].coeff2)*0.5;
			pair_coeff[i][j][0]=48 * epu*pow(sigma,12.0);
			pair_coeff[i][j][1]=24 * epu*pow(sigma,6.0);
		}
	}

}

void
Variables::ionRotation(void){
  random_device seed;
	double A,B,C,x,y,z;
	A=seed(),B=seed(),C=seed();
	for(auto &a : ions) {x=a.qx,y=a.qy,z=a.qz; ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z);}
}

void
Variables::ionInitialVelocity(double T) {
	random_device seed;
	default_random_engine engine(seed());
  for(auto &a : ions) {
    double matom=a.mass*1e-3/Nw;
    normal_distribution<> dist(0.0, sqrt(kb*T/matom));
		a.px=dist(engine)*1e-5;
		a.py=dist(engine)*1e-5;
		a.pz=dist(engine)*1e-5;
	}
}
