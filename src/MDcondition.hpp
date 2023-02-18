#pragma once
#include "constants.hpp"

//------------------------------------------------------------------------
class MDcondition {
public:

	MDcondition(void){
        MARGIN = 10.0;
        CUTOFF = 20.0;
        positionLogStep=0;
        logger=10000;
    };
	~MDcondition(void){};

    long int Noftimestep;
    long int step_relax;
    int positionLogStep;
    int logger;

    double L;
    double HL;
    double V;

    int Nof_around_gas;
    int Nof_around_vapor;

    double MARGIN;
    double CUTOFF;
    double CL2;
    double ML2;

private:
};
//------------------------------------------------------------------------
