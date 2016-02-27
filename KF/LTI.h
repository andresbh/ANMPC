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
    
    void append(
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> &A,
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> B)
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
    void append(
      Eigen::Matrix<_Scalar, Dynamic, 1> &A,
      Eigen::Matrix<_Scalar, Dynamic, 1> B)
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

    // createPredictionMatrices -----------------------
    void createPredictionMatrices(int Hp, int Hu) // prediction horizon, control horizon
    {
      this->Hp = Hp;
      this->Hu = Hu;
      Eigen::MatrixPower<Eigen::Matrix<_Scalar, Dynamic, Dynamic>> Apow(this->A);

      Eigen::Matrix<_Scalar, Dynamic, Dynamic> 
        sumApowB_Gamma = MatrixXd::Identity(this->D.rows(), this->D.cols());
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> 
        sumApowB_Theta = sumApowB_Gamma;
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> Gamma0, Theta0;
      Gamma0.resize(Hp*this->B.rows(), this->B.cols());   
      Theta0.resize(Hp*this->B.rows(), Hu*this->B.cols());
      Theta0.setZero();

      Phi.resize(Hp*this->C.rows(), this->C.cols());
      Gamma.resize(Hp*this->D.rows(), this->D.cols());
      Theta.resize(Hp*this->D.rows(), Hu*this->D.cols());

      Theta.setZero();
      
      int ii = 0, jj = 0;
      for (ii = 0; ii < Hp; ii++)
      { // populate Phi
        Phi.block(ii*this->C.rows(), 0, this->C.rows(), this->Phi.cols())
          = this->C*Apow(ii+1);//Apow1;
        
        if(ii==0)
        { // populate Gamma0
          Gamma0.block(ii*this->B.rows(), 0, this->B.rows(), Gamma0.cols())
            = this->B;
          
          // populate theta
          Theta0.block(ii*this->B.rows(), jj*this->B.cols(), this->B.rows(), this->B.cols())
            = this->B;
        }
        else
        { sumApowB_Gamma = //Apow2
          Apow(ii)*this->B;
          sumApowB_Theta = sumApowB_Gamma;
        
          // populate Gamma0
          Gamma0.block(ii*this->B.rows(), 0, this->B.rows(), Gamma0.cols())
            = (sumApowB_Gamma + Gamma0.block((ii - 1)*this->B.rows(), 0, this->B.rows(), Gamma0.cols()));
          
          // populate theta up to Hu
          Theta0.block(ii*this->B.rows(), 0, this->B.rows(), Gamma0.cols())
            << Gamma0.block(ii*this->B.rows(), 0, this->B.rows(), Gamma0.cols());

          Theta0.block(ii*this->B.rows(), this->B.cols(), this->B.rows(), Theta0.cols() - this->B.cols())
            << Theta0.block((ii-1)*this->B.rows(), 0, this->B.rows(), Theta0.cols() - this->B.cols());
          }
      }

      for (ii = 0; ii < Hp; ii++)
      {
        Gamma.block(ii*this->D.rows(), 0, this->D.rows(), this->Gamma.cols()) = 
          this->C*Gamma0.block(ii*this->B.rows(), 0, this->B.rows(), Gamma0.cols());
        for (jj = 0; jj < Hu; jj++)
        {
          Theta.block(ii*this->D.rows(), jj*this->D.cols(), this->D.rows(), this->D.cols())
            << this->C*Theta0.block(ii*this->B.rows(), jj*this->B.cols(), this->B.rows(), this->B.cols()); 
        }
      }
      
      //md.Th
      //  0.04332   0.00452   0.00000   0.00000   0.00000   0.00000
      //  0.00228   0.05422   0.00000   0.00000   0.00000   0.00000
      //  0.08416   0.01725   0.04332   0.00452   0.00000   0.00000
      //  0.00893   0.10699   0.00228   0.05422   0.00000   0.00000
      //  0.12268   0.03707   0.08416   0.01725   0.04332   0.00452
      //  0.01965   0.15835   0.00893   0.10699   0.00228   0.05422
      //  0.15900   0.06297   0.12268   0.03707   0.08416   0.01725
      //  0.03418   0.20833   0.01965   0.15835   0.00893   0.10699
      //  0.19326   0.09404   0.15900   0.06297   0.12268   0.03707
      //  0.05224   0.25698   0.03418   0.20833   0.01965   0.15835
      // std::cout << "theta:" << std::endl << this->Theta << std::endl;
      // std::cout << "gamma:" << std::endl << this->Gamma << std::endl;
      // std::cout << "phi:"   << std::endl << this->Phi   << std::endl;

    };

    // createOptimizationMatrices -----------------------
    void createOptimizationMatrices(
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> Q,
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> R,
      Eigen::Matrix<_Scalar, Dynamic, 1> umax,
      Eigen::Matrix<_Scalar, Dynamic, 1> umin,
      Eigen::Matrix<_Scalar, Dynamic, 1> dumax,
      Eigen::Matrix<_Scalar, Dynamic, 1> dumin
      )
    {
      Qq.resize(Theta.rows(),Theta.rows());
      Qq.setZero();
      Rr.resize(Theta.cols(), Theta.cols());
      Rr.setZero();

      Tc.resize(C.rows(),C.cols());
      Tc.setZero();
      Tc = this->C;

      for (int ii=0; ii < this->Hp; ii++)
      {
        Qq.block(ii*C.rows(), ii*C.rows(), C.rows(), C.rows()) << Q;
        this->blkdiag(Tc, C);
      }
      for (int ii = 0; ii < this->Hu; ii++)
      {
        Rr.block(ii*B.cols(), ii*B.cols(), B.cols(), B.cols()) << R;
      }
      H.resize(Rr.rows(),Rr.cols());
      H = Theta.transpose()*Qq*Theta + Rr;

      Eigen::Matrix<_Scalar, Dynamic, Dynamic> uf_block, uw_block, bu_block, au_block, z_block, qwef(2,1), qwe2(2, 1), qwew(2,1);
      uf_block.resize(2 * B.cols(), B.cols()); uf_block.setZero();
      uw_block.resize(2 * B.cols(), B.cols()); uw_block.setZero();
      bu_block.resize(2 * B.cols(), 1); bu_block.setZero();
      au_block.resize(2 * B.cols(), 1); au_block.setZero();
      z_block.resize(2 * B.cols(), 1); z_block.setZero(); // unused for zmin,zmax == -+inf
      for (int ii = 0; ii < B.cols(); ii++) // umin or umax == inf not possible
      {
        qwef(0, 0) = 1;qwef(1, 0) = -1;
        uf_block.block(B.cols() * ii, ii, B.cols(), 1) = qwef;
        qwe2(0, 0) = umax(ii); qwe2(1, 0) = -umin(ii);
        bu_block.block(B.cols() * ii, 0, B.cols(), 1) = qwe2;
        
        qwe2(0, 0) = dumax(ii); qwe2(1, 0) = -dumin(ii);
        au_block.block(B.cols() * ii, 0, B.cols(), 1) = qwe2;

        qwew(0, 0) = 1;qwew(1, 0) = -1;
        qwe2(0, 0) = dumax(ii); qwe2(1, 0) = -dumin(ii);
        uw_block.block(B.cols() * ii, ii, B.cols(), 1) = qwef;
      }
      
      F_con.resize((this->Hp+1)*this->B.cols(), this->Hu*this->B.cols() ); 
      F_con.setZero();
      
      f_con.resize(this->Hu*2*this->B.cols(), 1);
      f_con.setZero();

      for (int ii = 0; ii < this->Hu; ii++) // 
      { for (int jj = 0;jj <= ii;jj++)
        { F_con.block(uf_block.rows()*ii, jj*uf_block.cols(), uf_block.rows(), uf_block.cols()) =
            uf_block;         
        }
        f_con.block(bu_block.rows()*ii, 0, bu_block.rows(), bu_block.cols()) =
          bu_block;
      }

      W_con.resize((this->Hp + 1)*this->B.cols(), this->Hu*this->B.cols());
      W_con.setZero();
      w_con.resize(this->Hu * 2 * this->B.cols(), 1);
      w_con.setZero();

      for (int ii = 0; ii < this->Hu; ii++) // 
      {
        W_con.block(uw_block.rows()*ii, ii*uw_block.cols(), uw_block.rows(), uw_block.cols()) =
          uf_block;

        w_con.block(au_block.rows()*ii, 0, au_block.rows(), au_block.cols()) =
          au_block;
      }
    }

    // prediction --------------------------------------------
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> prediction
      (
        const Eigen::Matrix<_Scalar, Dynamic, Dynamic> u0, // prediction
        const Eigen::Matrix<_Scalar, Dynamic, Dynamic> Du0, // prediction
        const Eigen::Matrix<_Scalar, Dynamic, Dynamic> x0  // state
      )
    {
      Eigen::Matrix<_Scalar, Dynamic, Dynamic> Z0;
      Z0.resize(this->Gamma.rows(),1);

      Z0 = this->Phi*x0 +this->Gamma*u0.col(0) +this->Theta*Du0;

      //std::cout << "Z0:" << std::endl << Z0 << std::endl;

      return Z0;
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
    int Hp, Hu;
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> Phi; // prediction matrix
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> Gamma; // prediction matrix
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> Theta; // prediction matrix
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> Qq; // optimization matrix
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> Rr; // optimization matrix
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> H; // optimization matrix
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> G; // optimization matrix
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> Tc; // predi matrix
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> F_con; // optimization constraints
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> f_con; // optimization constraints
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> W_con; // optimization constraints
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> w_con; // optimization constraints
  };
} // namespace LTI_Object
#endif