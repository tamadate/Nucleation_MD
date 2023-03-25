# pragma once
#include "../../functions.hpp"


class Gyration : public Functions{
public:
    // check collision
	void postLoop(void);

    // initialization
    void initial(Variables *VARS, Physical *PP, MDcondition *CON, Observer *OBS);
    char filename[100];

	Gyration(int Step){
        type=0;
        step=Step;
    };
	~Gyration(void){};

private:

};
