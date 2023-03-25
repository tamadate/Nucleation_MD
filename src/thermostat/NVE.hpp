#pragma once
#include "../functions.hpp"

class NVE : public Functions {
	private:
		string statName="NVE";
	public:
		int type=1;
		void printName(void) {cout<<statName<<endl;}
		NVE(void){
			type=1;
			step=1;
		};
		~NVE(void){};
};

