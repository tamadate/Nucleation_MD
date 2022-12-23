#pragma once
#include "basePotential.hpp"
#include "AMBER/potentialAMBER.hpp"
#include "StillingerWeber/potentialSW.hpp"
#include "Tersoff/potentialTersoff.hpp"

class PotentialLJ : public Potential {
	private:
		string potName="LJ gas-ion";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialLJ(){};
		~PotentialLJ(){};
};

class PotentialLJCoulHybrid : public Potential {
	private:
		string potName="LJ vapor-vapor";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialLJCoulHybrid(){};
		~PotentialLJCoulHybrid(){};
};

class PotentialLJHybrid : public Potential {
	private:
		string potName="LJ vapor-gas";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialLJHybrid(){};
		~PotentialLJHybrid(){};
};

class PotentialLJCoul : public Potential {
	private:
		string potName="LJ vapor-ion";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialLJCoul(){};
		~PotentialLJCoul(){};
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
