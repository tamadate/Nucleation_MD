# pragma once
#include "variables.hpp"
#include "PhysicalProp.hpp"
#include "MDcondition.hpp"
#include "output/observer.hpp"

class Functions {
	private:
	public:
        Variables *vars;
        Physical *pp;
        MDcondition *con;
        Observer *obs;

        int type;
        int step;

        virtual void initial(Variables *VARS, Physical *PP, MDcondition *CON, Observer *OBS){
            vars=VARS;
            pp=PP;
            con=CON;
            obs=OBS;
        };

		virtual void preLoop(void) {}
		virtual void prePosition(void){};
        virtual void postPosition(void){};
        virtual void preForce(void){};
        virtual void postForce(void){};
        virtual void postLoop(void){};

        virtual void prePair(void){};
        virtual void postPair(void){};

		Functions(){};
		~Functions(){};
};
