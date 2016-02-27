#pragma once
#include <Eigen/Eigen>
#include "cppad/local/ad.hpp"
//#include <iostream>

using namespace Eigen;
using CppAD::AD;

class model //quadTank 
{friend class FG_info;
 public:
  int modelN_y;
  int modelN_u;
  bool integY_Ob;
  bool integU_Ob;
  
  model(void){}

  model(const int N_y0, const int N_z0, const int N_u0, const bool integY_b, const bool integU_b)
  { modelN_y=N_y0;
    modelN_u=N_u0;
    integY_Ob=integY_b;
    integU_Ob=integU_b;
    return;
  };

  int 
    modelcalc
      (const Matrix<AD<double>, Dynamic, 1> x, 
       const Matrix<AD<double>, Dynamic, 1> u, 
       const Matrix<AD<double>, Dynamic, 1> d, 
             Matrix<AD<double>, Dynamic, 1> &fvec)
  { double a1, a2, a3, a4;
    double A1, A2, A3, A4;
    double gamma1, gamma2;
    double k1, k2, g;
    a1 = 0.071; a2 = 0.071; a3 = 0.057; a4 = 0.057;
    A1 = 28.; A2 = 28.; A3 = 32.; A4 = 32.;
    gamma1 = 0.25; gamma2 = 0.35;
    k1 = 3.33; k2 = 3.35; g = 9.81;
    int i;

    //std::cout<<"xfunct: "<<x.transpose()<<std::endl;
    //std::cout<<"ufunct: "<<u.transpose()<<std::endl;
    //std::cout<<"dfunct: "<<d.transpose()<<std::endl;

    Matrix<AD<double>, Dynamic, 1> x_=x;
    for(i = 0; i<x_.rows(); i++ ) // sqrt non negative values
    { if(x_(i)<0)
        x_(i) = 0;
    }

    //std::cout<<"x_funct: "<<x_.transpose()<<std::endl;
    //std::cout<<"ufunct: "<<u.transpose()<<std::endl;
    //std::cout<<"dfunct: "<<d.transpose()<<std::endl;
    fvec(0) = -a1/A2*sqrt(2*g*x_(0)) +a3/A1*sqrt(2*g*x_(2)) +gamma1*k1/A1*(u(0)); 
    fvec(1) = -a2/A2*sqrt(2*g*x_(1)) +a4/A2*sqrt(2*g*x_(3)) +gamma2*k2/A2*(u(1)); 
    fvec(2) = -a3/A3*sqrt(2*g*x_(2)) +(1-gamma2)*k2/A3*u(1);
    fvec(3) = -a4/A4*sqrt(2*g*x_(3)) +(1-gamma1)*k1/A4*u(0);

    //std::cout<<"fv: "<<fvec.transpose()<<std::endl;
    return 1;
  }
  
};

///////////////////////////////////////////////////////////////////////////////