#include "../variables.hpp"

void ROTATION(double &X, double &Y, double &Z, double A, double B, double C, double x, double y, double z){
  X = cos(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+sin(C)*(y*cos(A)+z*sin(A));
  Y = -sin(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+cos(C)*(y*cos(A)+z*sin(A));
  Z = x*sin(B)-y*sin(A)*cos(B)+z*cos(A)*cos(B);
}
