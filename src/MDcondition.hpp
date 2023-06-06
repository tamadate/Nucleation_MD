#pragma once
#include "constants.hpp"

//------------------------------------------------------------------------
class MDcondition {
public:
    long int Noftimestep;   // number of iterations for main calculation
    long int step_relax;    // number of iterations for relaxation
    int positionLogStep;    // interval for recording vapor position
    int logger;             // interval for recording ion center properties
    int sampleStep;         // inverval for recording ""

    double L;   // simulation domain length
    double HL;  // simulation domain half length
    double V;   // simulation domain volume

    double dt;  // time step (delta t)

    int Nof_around_gas;     // number of gases in simulation domain
    int Nof_around_vapor;   // number of vapors in simulation domain

    double MARGIN;  // margin size for pair list
    double CUTOFF;  // cutoff length
    double CL2;     // square of cutoff
    double ML2;     // square of cutoff plus margin

    MDcondition(void){
        // setting default values
        MARGIN = 5.0;
        CUTOFF = 20.0;
        positionLogStep=0;
        logger=10000;
        sampleStep=10000;
        dt=0.5;
    };
	~MDcondition(void){};

private:
};
//------------------------------------------------------------------------
