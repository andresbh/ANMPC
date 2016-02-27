#pragma once
#ifndef COST_H
#define COST_H

//#include <omp.h>
#include <iostream>
#include <cppad/ipopt/solve.hpp>
#include "realModel.h"
#include "functorTank.h"
#include "EKF.h"

//#define NUMBER_THREADS 4

using CppAD::AD;

class COSTF: public model //quadTank
{
public:
  typedef AD<double> ADdouble;
  typedef CPPAD_TESTVECTOR( double ) Dvector;
  typedef CPPAD_TESTVECTOR( ADdouble ) ADvector;
  typedef Matrix< ADdouble , Dynamic, Dynamic> ADmatrix;

  // default constructor 
  COSTF();

  // default destructor
  virtual ~COSTF();

  COSTF(const int N_x0, const int N_y0, const int N_z0, const int N_u0, const int H_x0, const int H_u0, const VectorXi ublk0, const double h0);

  AD<double> costfct(const ADvector u);
  int system_output(Dvector u);
  int model_prediction(const Dvector u);
  virtual 
    int operator()(const Matrix<AD<double>, Dynamic, 1> x, 
                   const Matrix<AD<double>, Dynamic, 1> u, 
                   const Matrix<AD<double>, Dynamic, 1> d, 
                         Matrix<AD<double>, Dynamic, 1> &fvec); // virtual function from model class

  int N_x,  //number of states
      N_y,  //number of outputs
      N_z,  //number of controlled outputs
      N_u,  //number of inputs
      H_x,  //prediction horizon
      H_u;  //control horizon

  EKF<COSTF, AD<double>, AD<double>> *EKFS;

  Matrix<ADdouble, Dynamic, 1> measurements;
  Matrix<ADdouble, Dynamic, Dynamic> x_KF;
  Matrix<ADdouble, Dynamic, Dynamic> ref;
  Matrix<double, Dynamic, 1> u_applied;
  Matrix<ADdouble, Dynamic, Dynamic> Q;
  Matrix<ADdouble, Dynamic, Dynamic> R;
  Matrix<ADdouble, Dynamic, Dynamic> L;
  Matrix<int, 1, Dynamic> ublk;
    
};

#endif