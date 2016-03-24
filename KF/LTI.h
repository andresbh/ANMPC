#pragma once
#ifndef LTI_H
#define LTI_H
#include <iostream>
#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>

namespace LTI_Object {
  using namespace Eigen;
  template <typename _Scalar>
  class LTI
  {protected:
    template<typename T>
    Eigen::Matrix<T, Dynamic, 1> cumsum(const Eigen::Matrix<T, Dynamic, 1> vec) {
      Eigen::Matrix<T, Dynamic, 1> temp;

      temp.resizeLike(vec);
      temp = vec;
      for (int j = 1; j < vec.rows(); ++j)
        temp(j, 0) += temp(j - 1, 0);
      return temp;
    };

    void blkdiag(
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> &A,
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> B
      )
    {
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> Temp;
      Temp.resize(A.rows()+B.rows(), A.cols() + B.cols());
      Temp.setZero();
      Temp.block(0, 0, A.rows(), A.cols()) = A;
      Temp.block(A.rows(), A.cols(), B.rows(), B.cols()) = B;
      A.resize(Temp.rows(),Temp.cols());
      A = Temp;
    }

    Eigen::Matrix<_Scalar, Dynamic, Dynamic> appendvr( //vertical append
      const Eigen::Matrix<_Scalar, Dynamic, Dynamic> A,
      const Eigen::Matrix<_Scalar, Dynamic, Dynamic> B)
    {
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> temp;
      temp.resizeLike(A);
      temp = A;
      appendv(temp, B);
      return temp;
    };

    void appendv( //vertical append
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> &A,
      const Eigen::Matrix<_Scalar, Dynamic, Dynamic> B)
    {
      if (B.cols() != A.cols())
      {
        std::cout << "Append function error matrix sizes. Size A: " << A.cols() << " size B: " << B.cols() << std::endl;
        assert(B.cols() != A.cols());
      }
        
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> Temp;
      Temp.resize(A.rows() + B.rows(), A.cols());
      Temp.setZero();
      Temp.block(0, 0, A.rows(), A.cols()) = A;
      Temp.block(A.rows(), 0, B.rows(), B.cols()) = B;
      A.resize(Temp.rows(), Temp.cols());
      A = Temp;
    };
    
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> appendhr( //horizontal append
      const Eigen::Matrix<_Scalar, Dynamic, Dynamic> A,
      const Eigen::Matrix<_Scalar, Dynamic, Dynamic> B)
    {
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> temp;
      temp.resizeLike(A);
      temp = A;
      appendh(temp, B);
      return temp;
    };

    void appendh( //horizontal append
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> &A,
      const Eigen::Matrix<_Scalar, Dynamic, Dynamic> B)
    {
      if (B.rows() != A.rows())
      {
        std::cout << "Append function error matrix sizes. Size A: " << A.rows() << " size B: " << B.rows() << std::endl;
        assert(B.rows() != A.rows());
      }

      Eigen::Matrix<_Scalar, Dynamic, Dynamic> Temp;
      Temp.resize(A.rows(), A.cols()+B.cols());
      Temp.setZero();
      Temp.block(0, 0, A.rows(), A.cols()) = A;
      Temp.block(0, A.cols(), B.rows(), B.cols()) = B;
      A.resize(Temp.rows(), Temp.cols());
      A = Temp;
    };

    template<typename T>
    void append(
      Eigen::Matrix<T, Dynamic, 1> &A,
      const Eigen::Matrix<T, Dynamic, 1> B
      )
    {
      if (B.cols() != A.cols())
      {
        std::cout << "Append function error vector sizes. Size A: " << A.cols() << " size B: " << B.cols() << std::endl;
        assert(B.cols() != A.cols());
      }

      Eigen::Matrix<T, Dynamic, 1> Temp;
      Temp.resize(A.rows() + B.rows(), A.cols());
      Temp.setZero();
      Temp.block(0, 0, A.rows(), A.cols()) = A;
      Temp.block(A.rows(), 0, B.rows(), B.cols()) = B;
      A.resize(Temp.rows(), Temp.cols());
      A = Temp;

    };

   public:
    //controllable canonical form TF
    LTI(void)
    {
    };
    ~LTI(void)
    { 
    };

    LTI(
      Eigen::Matrix<_Scalar, Dynamic, 1> y0,
      Eigen::Matrix<_Scalar, 1, Dynamic> A0,
      Eigen::Matrix<_Scalar, 1, Dynamic> B0
      )
    { if (A0.size() < B0.size())
        return;
      A.setZero(A0.size(), A0.size());
      A.topRightCorner(A0.size() - 1, A0.size() - 1).setIdentity(A0.size() - 1, A0.size() - 1);
      A.row(A.rows() - 1) = -A0;

      B.setZero(A0.size(), 1);
      B(A0.size() - 1, 0) = 1;

      C.setZero(1, A0.size());
      C.bottomLeftCorner(1, B0.size()) = B0;

      D.setZero(1, 1);

      x.setZero(A.rows());
      y = y0;
      //std::cout<<"A: "<<A<<std::endl
      //  <<"B: "<<B<<std::endl
      //  <<"C: "<<C<<std::endl
      //  <<"D: "<<D<<std::endl
      //  <<"x: "<<x<<std::endl;
    };

    // definition of class ---------------------
    LTI(
      Eigen::Matrix<_Scalar, Dynamic, 1>       y0,
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> A0,
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> B0,
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> C0,
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> D0
      )
    { //if(A0.size()<=B0.size())
      //  return;
      A = A0;
      B = B0;
      C = C0;
      D = D0;

      x.setZero(A.rows());
    };

    virtual Eigen::Matrix<_Scalar, Dynamic, 1> update(const Eigen::Matrix<_Scalar, Dynamic, 1> u0)
    { x = A*x + B*u0;
      return C*x + D*u0;
    };

        
    friend LTI<_Scalar> operator*(
       const LTI<_Scalar> SYS2
      )
    { return this->series(SYS2);
    };

    LTI<_Scalar> series(
      const LTI<_Scalar> SYS2
      )
    { LTI<_Scalar> LTItemp;
      Matrix<_Scalar, Dynamic, Dynamic> MatrixTemp;
      MatrixTemp = this->B*this->C;
      LTItemp.A = MatrixXd::Zero(this->A.rows() + SYS2.A.rows(), this->A.cols() + SYS2.A.cols());
      LTItemp.A.topLeftCorner(this->A.rows(), this->A.cols()) = this->A;
      LTItemp.A.bottomRightCorner(this->A.rows(), this->A.cols()) = SYS2.A;
      LTItemp.A.bottomLeftCorner(MatrixTemp.rows(), MatrixTemp.cols()) = MatrixTemp;

      MatrixTemp = this->D*SYS2.B;
      LTItemp.B = MatrixXd::Zero(this->B.rows() + MatrixTemp.rows(), this->B.cols());
      LTItemp.B.block(0, 0, this->B.rows(), this->B.cols()) = this->B;
      LTItemp.B.block(this->B.rows() + 1, 0, MatrixTemp.rows(), this->B.cols()) = MatrixTemp;

      MatrixTemp = SYS2.D*this->C;
      LTItemp.C = MatrixXd::Zero(SYS2.C.rows(), MatrixTemp.cols() + SYS2.C.cols());
      LTItemp.C.block(0, 0, SYS2.C.rows(), MatrixTemp.rows()) = MatrixTemp;
      LTItemp.C.block(0, MatrixTemp.rows() + 1, SYS2.C.rows(), this->B.rows()) = SYS2.C;

      MatrixTemp = this->D*SYS2.D;
      LTItemp.D = MatrixTemp;

      return LTItemp;
    };

    friend LTI<_Scalar> operator+(
      const LTI<_Scalar> SYS2
      )
    { return this->parallel(SYS2);
    };

    LTI<_Scalar> parallel(
      const LTI<_Scalar> SYS2
      )
    {
      LTI<_Scalar> temp;
      return temp;
    };

    LTI<_Scalar> feedback(
      const LTI<_Scalar> SYS1,
      const LTI<_Scalar> SYS2
      )
    { LTI<_Scalar> temp;
      return temp;
    };

    Eigen::Matrix<_Scalar, Dynamic, Dynamic> A;
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> B;
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> C;
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> D;
    Eigen::Matrix<_Scalar, Dynamic, 1>       y;
    Eigen::Matrix<_Scalar, Dynamic, 1>       x;

  };
} // namespace LTI_Object
#endif