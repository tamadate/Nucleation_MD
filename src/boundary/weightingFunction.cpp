#include "../md.hpp"

double
MD:: weightFunc(double dr2){
  double returnValue=0;
  if(dr2<RCG2) {
    returnValue=1;
  } else if(dr2>RAA2) {
    returnValue=0;
  } else {
    double dr=(sqrt(dr2)-RCG);
    double C=cos(M_PI/BoundL*0.25*dr);
    returnValue=C*C;
  }
  return returnValue;

}
