#pragma once
#ifndef LKF_H
#define LKF_H
#include "LTI.h"
//#include "../NMPCipoptccpAD/Cost.hpp"
#include <iostream>

using namespace Eigen;

template <typename _Scalar>
class LKF: public LTI< _Scalar > 
{public:
  LKF(void){};
  
  // state space form
  LKF(
     Matrix<_Scalar, Dynamic, 1>       x0,
     Matrix<_Scalar, Dynamic, 1>       y0,
     Matrix<_Scalar, Dynamic, Dynamic> P0,
     Matrix<_Scalar, Dynamic, Dynamic> Q0,
     Matrix<_Scalar, Dynamic, Dynamic> R0,
     Matrix<_Scalar,       1, Dynamic> A0,
     Matrix<_Scalar,       1, Dynamic> B0
    ): LKF::LTI(y0, A0, B0)
{ x = x0;
  K.resize( C.rows(), A.rows());
  K.setZero( C.rows(), A.rows());
  P = P0;
  Q = Q0;
  R = R0;
  y.resize( C.rows(), 1);
  return;
};

  LKF(
     Matrix<_Scalar, Dynamic, 1>       x0,
     Matrix<_Scalar, Dynamic, 1>       y0,
     Matrix<_Scalar, Dynamic, Dynamic> P0,
     Matrix<_Scalar, Dynamic, Dynamic> Q0,
     Matrix<_Scalar, Dynamic, Dynamic> R0,
     Matrix<_Scalar, Dynamic, Dynamic> A0,
     Matrix<_Scalar, Dynamic, Dynamic> B0,
     Matrix<_Scalar, Dynamic, Dynamic> C0,
     Matrix<_Scalar, Dynamic, Dynamic> D0
    ) : LKF::LTI(y0, A0, B0, C0, D0)
{ x = x0;
  K.resize( C.rows(), A.rows());
  K.setZero( C.rows(), A.rows());
  P = P0;
  Q = Q0;
  R = R0;
  y.resize( C.rows(), 1);
  return;
};

  Matrix<_Scalar, Dynamic, 1> update(const Matrix<_Scalar, Dynamic, 1> y0, const Matrix<_Scalar, Dynamic, 1> u0)
{ //std::cout<<"\nx: "<< x;
  x =A*x +B*u0;
  y =C*x +D*u0;

  //std::cout<<"\ny: "<< y<<std::endl;
  //std::cout<<"\nx: "<< x<<std::endl;
  P = A*P*A.transpose() + Q;
  K = ((C*P*C.transpose() + R).transpose()).lu().solve((P*C.transpose()).transpose()).transpose();
  //std::cout<<"\n  e: "<< (y-C*x)<< "\t y: "<<y <<"\t Cx: "<<C*x <<"\t u: "<<u.transpose()<<std::endl;
  x+= K*(y0-y);
  P+= -K*C*P;
  return x;
};

  // vars
  Matrix<_Scalar, Dynamic, Dynamic> P;
  Matrix<_Scalar, Dynamic, Dynamic> Q;
  Matrix<_Scalar, Dynamic, Dynamic> R;

protected:
  // methods

  // vars
  Matrix<_Scalar, Dynamic, Dynamic> K;
};

#endif