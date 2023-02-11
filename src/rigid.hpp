#pragma once
#include "variables.hpp"
#include "array"

//------------------------------------------------------------------------
class Rigid {
public:

	Rigid(void);
	~Rigid(void){

    };

    std::vector<Bond> rigid_pairs_ion;
    std::vector<vector<Bond>> rigid_pairs_vapor;
    std::vector<vector<Bond>> rigid_pairs_gas;
    std::vector<std::array<double,3>> r;

    int MaxIteration;

    void recordOldPosition(std::vector<Atom> &target);
    void updatePosition(Variables *vars,double dt);
    void shake_ion(Variables *vars);
    void shake_vapor(Variables *vars,int vaporID);
    void shake_gas(Variables *vars,int vaporID);


private:
};
//------------------------------------------------------------------------
