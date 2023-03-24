#include "../variables.hpp"


class Collision{
public:
    Variables *vars;
    std::vector<long int> collisionFlagVapor;
	void checkCollision(int itime);
    void outputOut(int i);
    void outputIn(int i);
    void init(Variables *VARS){
        vars=VARS;
        int Nvapor=vars->vapors.size();
        collisionFlagVapor.resize(Nvapor);
        for (auto i : collisionFlagVapor) i=0;
    };

	Collision(){};
	~Collision(void){};

private:

};
