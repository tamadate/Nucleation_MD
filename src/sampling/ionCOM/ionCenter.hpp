# pragma once
#include "../../functions.hpp"


class IonCenter : public Functions{
public:
    // check collision
	void postLoop(void);

    // initialization
    void initial(Variables *VARS, Physical *PP, MDcondition *CON, Observer *OBS);
    char filename[100];

	IonCenter(int Step){
        type=0;
        step=Step;
    };
	~IonCenter(void){};

private:

};
