#include "realModel.h"
#include <iostream>

int realModel(double *fvec, const double *u, const double* x)
{ double a1 = 0.071, a2 = 0.071, a3 = 0.057, a4 = 0.057;
  double A1 = 28, A2 = 28, A3 = 32, A4 = 32;
  double gamma1 = 0.25, gamma2 = 0.35;
  double k1 = 3.33, k2 = 3.35;
  double g = 9.81;
  for(int i = 0; i<4; i++ ) // sqrt non negative values
  { if(x[i]<0)
      fvec[i] = 0;
    else
      fvec[i] = x[i];
  }
  //std::cout<<"u_1: "<<u[0]<<" u_2: "<<u[1]<<std::endl;

  fvec[0] = -a1/A2*sqrt(2*g*fvec[0]) +a3/A1*sqrt(2*g*fvec[2]) +gamma1*k1/A1*u[0];
  fvec[1] = -a2/A2*sqrt(2*g*fvec[1]) +a4/A2*sqrt(2*g*fvec[3]) +gamma2*k2/A2*u[1];
  fvec[2] = -a3/A3*sqrt(2*g*fvec[2]) +(1-gamma2)*k2/A3*u[1];
  fvec[3] = -a4/A4*sqrt(2*g*fvec[3]) +(1-gamma1)*k1/A4*u[0];
  
  return 0;
}

