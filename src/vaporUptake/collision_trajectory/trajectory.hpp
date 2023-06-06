# pragma once
#include "../../functions.hpp"


class Trajectory : public Functions{
public:
    // check collision
	void postPosition(void);

    // initialization
    void initial(Variables *VARS, Physical *PP, MDcondition *CON, Observer *OBS);

    char filename[100];
    double da;  // arrival (collision) distance
    double dl;  // departure (leaving) distance
    std::vector<vector<double>> sigma216;  // collision distance (2^(1/6)*sigma)
    std::vector<long int> collisionFlagVapor;

	Trajectory(double arrival, double leave, int Step) {
        da=arrival;
        dl=leave;
        type=0;
        step=Step;
    };
	~Trajectory(void){};

private:

};
