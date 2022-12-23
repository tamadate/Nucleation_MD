#include "../variables.hpp"

Variables::Variables(void) {
  time = 0.0;
  Utotal.Uion=Utotal.Ugas=Utotal.Uvap=Utotal.Ugi=Utotal.Ugg=Utotal.Uvg=Utotal.Uvi=Utotal.Uvv=0;

  Molecule mol;
  mol.qx=mol.qy=mol.qz;
  mol.px=mol.py=mol.pz;
  mol.fx=mol.fy=mol.fz;
  MolID.resize(3);
  MolID[0].push_back(0);
}


void
Variables::read_initial(char* ionFile, char* vaporFile) {
  int Ntype_ion=readIonFile(ionFile);
  int Ntype_vapor=readVaporFile(vaporFile);
  setBMHPotential();
  setCrossPotentials(Ntype_ion,Ntype_vapor);
  setInitialRegion();
  ionRotation();
}


void
Variables::setCrossPotentials(int Nion,int Nvapor){
  int Natypes=atypes.size();
  pair_coeff.resize(Natypes);
  for (int i=0;i<Natypes;i++){
    pair_coeff[i].resize(Natypes);
    for (int j=0;j<Natypes;j++){
      pair_coeff[i][j].resize(2);
    }
  }
  for (int i=0;i<Natypes;i++){
    for (int j=0;j<Natypes;j++){
      double epu=sqrt(atypes[i].coeff1*atypes[j].coeff1);
      double sigma=(atypes[i].coeff2+atypes[j].coeff2)*0.5;
      //sigma=sqrt(atypes_v[i].coeff2*atypes_v[j].coeff2);
      pair_coeff[i][j][0]=48 * epu*pow(sigma,12.0);
      pair_coeff[i][j][1]=24 * epu*pow(sigma,6.0);
    }
  }

  double epu11=0.14397; // gas CG potential
  double sigma11=3.798; // gas CG potential
  double epu22=0.14397; // vapor CG potential
  double sigma22=3.798; // vapor CG potential
  double epu12=sqrt(epu11*epu22);
  double sigma12=(sigma11+sigma11)*0.5;
  pair_coeff_CG[1][1][0]= 48 * epu11*pow(sigma11,12.0);
  pair_coeff_CG[1][1][1]= 24 * epu11*pow(sigma11,6.0);
  pair_coeff_CG[2][2][0]= 48 * epu22*pow(sigma22,12.0);
  pair_coeff_CG[2][2][1]= 24 * epu22*pow(sigma22,6.0);
  pair_coeff_CG[1][2][0]=pair_coeff_CG[2][1][1]= 48 * epu12*pow(sigma12,12.0);
  pair_coeff_CG[1][2][1]=pair_coeff_CG[2][1][1]= 24 * epu12*pow(sigma12,6.0);
}

void
Variables::setBMHPotential(void){
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


void
Variables::ionRotation(void){
  random_device seed;
	double A,B,C,x,y,z;
	A=seed(),B=seed(),C=seed();
	for(auto &a : Molecules[0].inAtoms) {x=a.qx,y=a.qy,z=a.qz; ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z);}
    int is=Molecules[0].inAtoms.size();
    for(int i=0;i<is-1;i++) {
        for(int j=i+1;j<is;j++){
            int flag=0;
            for (auto &d : dihedrals) {
                int I=d.atom1;
                int J=d.atom2;
                int K=d.atom3;
                int L=d.atom4;
                if(i==I){if(j==J||j==K||j==L) flag=1;}
                if(i==J){if(j==I||j==K||j==L) flag=1;}
                if(i==K){if(j==J||j==I||j==L) flag=1;}
                if(i==L){if(j==J||j==K||j==I) flag=1;}
            }
            if (flag==0){
                Pair p;
                p.i=i;
                p.j=j;
                ion_pairs.push_back(p);
            }
        }
    }
}

void
Variables::ionInitialVelocity(double T) {
	random_device seed;
	default_random_engine engine(seed());
  for(auto &a : Molecules[0].inAtoms) {
    double matom=a.mass*1e-3/Nw;
    normal_distribution<> dist(0.0, sqrt(kb*T/matom));
		a.px=dist(engine)*1e-5;
		a.py=dist(engine)*1e-5;
		a.pz=dist(engine)*1e-5;
	}
}

void
Variables::setInitialRegion(void) {
  int Nmol=Molecules.size();
  Region.resize(Nmol);
  int i=0;
  for (auto &mol : Molecules){
    double dx=mol.qx-Molecules[0].qx;
    double dy=mol.qy-Molecules[0].qy;
    double dz=mol.qz-Molecules[0].qz;
    double dr2=dx*dx+dy*dy+dz*dz;
    if(dr2<RO2) Region[i]=00000001;
    else {
      if(dr2<RI2) Region[i]=00000011;
      else Region[i]=00000010;
    }
    i++;
  }
}
