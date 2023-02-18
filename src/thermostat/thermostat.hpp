#pragma once
#include "../variables.hpp"
#include "../output/observer.hpp"

//------------------------------------------------------------------------
class Thermostat {
	private:
		string statName="NVE";
	public:
		virtual void printName(void) {cout<<statName<<endl;}
		virtual void tempControl(double dt){};
		virtual void computeZeta(double dt){};
		virtual void initial(){};
		Thermostat(){};
		~Thermostat(){};
};

class ThermostatNH : public Thermostat {
	private:
		string statName="NVT Nose-Hoover";
	public:
        Variables *vars;
        Observer *obs;
        double zeta;
        double Ttarget;
        double Q_inv;
		void printName(void) {cout<<statName<<endl;}
		void tempControl(double dt);
		void computeZeta(double dt);
		void initial(){};
		ThermostatNH(Variables *VARS, Observer *OBS, double t){
            vars=VARS;
            obs=OBS;
            Ttarget=t;
            Q_inv = 0.0001;
        };
		~ThermostatNH(){};
};

class ThermostatVscale : public Thermostat {
	private:
		string statName="NTV Velocity scaling";
	public:
        Variables *vars;
        Observer *obs;
        double zeta;
        double Ttarget;
		void printName(void) {cout<<statName<<endl;}
		void tempControl(double dt);
		void computeZeta(double dt){};
		void initial(){};
		ThermostatVscale(Variables *VARS, Observer *OBS, double t){
            vars=VARS;
            obs=OBS;
            Ttarget=t;
        };
		~ThermostatVscale(){};
};