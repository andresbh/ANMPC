#pragma once

#include <Eigen/Eigen>

using namespace Eigen;
///////////////////////////////////////////////////////////////////////////////
class model
{public:
 model(void){};
 int operator()(const Matrix<double, Dynamic, 1> x, 
                const Matrix<double, Dynamic, 1> u, 
                const Matrix<double, Dynamic, 1> d,
                      Matrix<double, Dynamic, 1> &fvec, 
                      double m=0.25, double M=1, double l=1.0, double theta = 0.0)
 { 
   // inverse pendulum
   double g=9.81; 
   fvec(0) = x(1);  //theta'
   fvec(1) = cos(x(0)+theta)*( g*m*cos(x(0)+theta)*sin(x(0)+theta)-m*l*x(1)*x(1)*sin(x(0)+theta) +u(0))/(l*(M+m*sin(x(0)+theta)*sin(x(0)+theta))) +g/l*sin(x(0)+theta) -0.7*x(1); //theta''
   fvec(2) = x(3);  //x'
   fvec(3) = (g*m*cos(x(0)+theta)*sin(x(0)+theta)-m*l*x(1)*x(1)*sin(x(0)+theta) +u(0))/(M+m*sin(x(0)+theta)*sin(x(0)+theta));//x''

   return 0;
 };

 void processDynamics(const Matrix<double, Dynamic, 1> x, 
                      const Matrix<double, Dynamic, 1> u, 
                      const Matrix<double, Dynamic, 1> d,
                            Matrix<double, Dynamic, 1> &fvec)
 { this->operator()(x, u, d, fvec);
 };

 void processOutput(const Matrix<double, Dynamic, 1> x, 
                    const Matrix<double, Dynamic, 1> u, 
                    const Matrix<double, Dynamic, 1> d,
                          Matrix<double, Dynamic, 1> &fvec)
 { fvec(0) = x(2);
   fvec(1) = x(3);
 };
};
///////////////////////////////////////////////////////////////////////////////