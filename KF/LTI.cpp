#include "LTI.h"
//#include <iostream>

template <typename _Scalar> 
LTI<_Scalar>::LTI(void)
{
}

template <typename _Scalar> 
LTI<_Scalar>::~LTI(void)
{
}

template <typename _Scalar> 
LTI<_Scalar>::LTI(const Matrix<_Scalar, Dynamic, 1>       y0,
                  const Matrix<double ,       1, Dynamic> A0,
                  const Matrix<double ,       1, Dynamic> B0
 )
{ if(A0.size()<B0.size())
    return;
  A.setZero(A0.size(),A0.size());
  A.topRightCorner( A0.size()-1, A0.size()-1).setIdentity( A0.size()-1, A0.size()-1);
  A.row(A.rows()-1) = -A0;
  
  B.setZero(A0.size(),1);
  B(A0.size()-1,0) = 1;

  C.setZero(1,A0.size());
  C.bottomLeftCorner(1, B0.size()) = B0;

  D.setZero(1,1);

  x.setZero(A.rows());
  y = y0;
  //std::cout<<"A: "<<A<<std::endl
  //  <<"B: "<<B<<std::endl
  //  <<"C: "<<C<<std::endl
  //  <<"D: "<<D<<std::endl
  //  <<"x: "<<x<<std::endl;
  
}

template <typename _Scalar> 
LTI<_Scalar>::LTI(const Matrix<_Scalar, Dynamic, 1>       y0,
                  const Matrix<double , Dynamic, Dynamic> A0,
                  const Matrix<double , Dynamic, Dynamic> B0,
                  const Matrix<double , Dynamic, Dynamic> C0,
                  const Matrix<double , Dynamic, Dynamic> D0
 )
{ //if(A0.size()<=B0.size())
  //  return;
  A=A0;
  B=B0;
  C=C0;
  D=D0;

  x.setZero(A.rows());
}

template <typename _Scalar> 
Matrix<_Scalar, Dynamic, 1> LTI<_Scalar>::update(const Matrix<_Scalar, Dynamic, 1> u0)
{ x = A*x +B*u0;
  return C*x +D*u0;
}

template <typename _Scalar> 
LTI<_Scalar> operator*(    
    const LTI<_Scalar> SYS1,
    const LTI<_Scalar> SYS2
)
{ LTI<_Scalar> temp;
  return temp.series(SYS1, SYS2);
}

template <typename _Scalar> 
LTI<_Scalar> LTI<_Scalar>::series(    
    const LTI<_Scalar> SYS1,
    const LTI<_Scalar> SYS2
)
{ LTI<_Scalar> LTItemp;
  Matrix<_Scalar, Dynamic, Dynamic> MatrixTemp;
  MatrixTemp = SYS1.B*SYS1.C;
  LTItemp.A = MatrixXd::Zero(SYS1.A.rows()+SYS2.A.rows(), SYS1.A.cols()+SYS2.A.cols()); 
  LTItemp.A.topLeftCorner(SYS1.A.rows(),SYS1.A.cols()) = SYS1.A;
  LTItemp.A.bottomRightCorner(SYS1.A.rows(),SYS1.A.cols()) = SYS2.A;
  LTItemp.A.bottomLeftCorner(MatrixTemp.rows(),MatrixTemp.cols()) = MatrixTemp;
  
  MatrixTemp = SYS1.D*SYS2.B;
  LTItemp.B = MatrixXd::Zero(SYS1.B.rows()+MatrixTemp.rows(), SYS1.B.cols()); 
  LTItemp.B.block(0, 0, SYS1.B.rows(), SYS1.B.cols()) = SYS1.B;
  LTItemp.B.block(SYS1.B.rows()+1, 0, MatrixTemp.rows(), SYS1.B.cols()) = MatrixTemp;

  MatrixTemp = SYS2.D*SYS1.C;
  LTItemp.C = MatrixXd::Zero(SYS2.C.rows(), MatrixTemp.cols()+SYS2.C.cols()); 
  LTItemp.C.block(0, 0, SYS2.C.rows(), MatrixTemp.rows()) = MatrixTemp;
  LTItemp.C.block(0, MatrixTemp.rows()+1, SYS2.C.rows(), SYS1.B.rows()) = SYS2.C;

  MatrixTemp = SYS1.D*SYS2.D;
  LTItemp.D = MatrixTemp;

  return LTItemp;
}

template <typename _Scalar> 
LTI<_Scalar> operator+(    
    const LTI<_Scalar> SYS1,
    const LTI<_Scalar> SYS2
)
{ LTI<_Scalar> temp;
  return temp.parallel(SYS1, SYS2);
}

template <typename _Scalar> 
LTI<_Scalar> LTI<_Scalar>::parallel(    
    const LTI<_Scalar> SYS1,
    const LTI<_Scalar> SYS2
)
{ LTI<_Scalar> temp;
  return temp;
}

template <typename _Scalar> 
LTI<_Scalar> LTI<_Scalar>::feedback(
    const LTI<_Scalar> SYS1,
    const LTI<_Scalar> SYS2
)
{ LTI<_Scalar> temp;
  return temp;
}

//The explicit instantiation part
template class LTI< double >; 
template class LKF< COSTF::ADdouble >;