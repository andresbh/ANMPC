#include "LKF.h"

template <typename _Scalar> 
LKF<_Scalar>::LKF(void)
{}

template <typename _Scalar> 
LKF<_Scalar>::LKF(
    const Matrix<_Scalar, Dynamic, 1>       x0,
    const Matrix<_Scalar, Dynamic, 1>       y0,
    const Matrix<double , Dynamic, Dynamic> P0,
    const Matrix<double , Dynamic, Dynamic> Q0,
    const Matrix<double , Dynamic, Dynamic> R0,
    const Matrix<double ,       1, Dynamic> A0,
    const Matrix<double ,       1, Dynamic> B0
    ) : LKF::LTI(y0, A0, B0)
{ x = x0;
  K.resize( C.rows(), A.rows());
  K.setZero( C.rows(), A.rows());
  P = P0;
  Q = Q0;
  R = R0;
  y.resize( C.rows(), 1);
  return;
}

template <typename _Scalar> 
LKF<_Scalar>::LKF(
    const Matrix<_Scalar, Dynamic, 1>       x0,
    const Matrix<_Scalar, Dynamic, 1>       y0,
    const Matrix<double , Dynamic, Dynamic> P0,
    const Matrix<double , Dynamic, Dynamic> Q0,
    const Matrix<double , Dynamic, Dynamic> R0,
    const Matrix<double , Dynamic, Dynamic> A0,
    const Matrix<double , Dynamic, Dynamic> B0,
    const Matrix<double , Dynamic, Dynamic> C0,
    const Matrix<double , Dynamic, Dynamic> D0
   ) : LKF::LTI(y0, A0, B0, C0, D0)
{ x = x0;
  K.resize( C.rows(), A.rows());
  K.setZero( C.rows(), A.rows());
  P = P0;
  Q = Q0;
  R = R0;
  y.resize( C.rows(), 1);
  return;
}

template <typename _Scalar> 
Matrix<_Scalar, Dynamic, 1> LKF<_Scalar>::update(const Matrix<_Scalar, Dynamic, 1> y0, const Matrix<_Scalar, Dynamic, 1> u0)
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
}

//The explicit instantiation part
template class LKF<double>;
template class LKF< AD<double> >;
