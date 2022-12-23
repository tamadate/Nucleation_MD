#pragma once
#include "../variables.hpp"

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
