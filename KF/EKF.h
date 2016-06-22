#pragma once
#ifndef EKF_H
#define EKF_H
#include "LKF.h"
#include <iostream>
//#include "../NMPCipoptccpAD/Cost.hpp"

template <class _Functor, typename _Scalar, typename _EKFScalar >
class EKF: public LKF< _Scalar >, public _Functor
{ public:
   EKF( //const _Functor& F0,
         const Matrix<_Scalar, Dynamic, 1> x0,
         const Matrix<_Scalar, Dynamic, 1> y0,
         const Matrix<_Scalar, Dynamic, 1> u0,
         const double h0,
         const Matrix<_Scalar, Dynamic, Dynamic> C0, // outputs
         const int N_x0, // real number of states
         const int N_u0,
         const int N_z0, // controlled outputs
         const Matrix<_Scalar, Dynamic, Dynamic> P0, // cov matrix estimation
         const Matrix<_Scalar, Dynamic, Dynamic> Q0, // weight for outputs
         const Matrix<_Scalar, Dynamic, Dynamic> R0  // weight for inputs
      ): EKF::LKF( x0, y0, P0, Q0, R0, 
                   MatrixXd::Zero(C0.cols(),C0.cols()).cast<_EKFScalar>(), //A AD<double>
                   MatrixXd::Zero(C0.cols(),u0.rows()).cast<_EKFScalar>(), //B
                   C0.cast<_EKFScalar>(),                                  //C
                   MatrixXd::Zero(C0.rows(),u0.rows()).cast<_EKFScalar>()  //D
                 ),
         N_x(N_x0),
         N_u(N_u0),
         N_z(N_z0)
  {   h = h0;
      x = x0;
      dx0_.resize(x.size());
      x0_.resize(x.size());
  };
   
   int operator()( const Matrix<_Scalar, Dynamic, 1> x0 )
   {  x = x0;
      return 0;
   };

   void update( const Matrix<_Scalar, Dynamic, 1> y0, 
                const Matrix<_Scalar, Dynamic, 1> u0
              )
   { //std::cout<<"xkf_in: "<<x.transpose()<<std::endl;
     updateJ(A, y0, u0, true); // continuous state space
     A = MatrixXd::Identity(A.rows(), A.cols()).cast<_EKFScalar>() +h*A;
     //std::cout<<"A_kf: "<<std::endl<<A<<std::endl;
     updateJ(B, y0, u0, false); // continuous state space
     B = h*B;
     //std::cout<<"B_kf: "<<std::endl<<B<<std::endl;
     x0_ = x;

     LKF::update(y0, u0);

     _Functor::operator()(x_.head(N_x),
                          u_,
                          x_.segment(N_x, N_z),
                          dx0_);

     x = (x_ +dx0_*h).cast<_Scalar>();
     //std::cout<<"xkf_out: "<<x.transpose()<<std::endl;
     return;
   };

   Matrix<_Scalar, Dynamic, 1> 
    update(       Matrix<_Scalar, Dynamic, 1> xkf,
                  const Matrix<_Scalar, Dynamic, 1> y0,
                  const Matrix<_Scalar, Dynamic, 1> u0
              )
   { //std::cout<<"xkf_in: "<<xkf.transpose()<<std::endl;
     updateJ(A, y0, u0, true); // continuous state space
     A = MatrixXd::Identity(A.rows(), A.cols()).cast<_EKFScalar>() +h*A;
     //std::cout<<"A_kf: "<<std::endl<<A<<std::endl;
     updateJ(B, y0, u0, false); // continuous state space
     B = h*B;
     //std::cout<<"B_kf: "<<std::endl<<B<<std::endl;
     x0_ = xkf;

     LKF::update(y0, u0);

     _Functor::operator()(x_.head(N_x),
                          u_,
                          x_.segment(N_x, N_z),
                          dx0_);

     return (x_ +dx0_*h).cast<_Scalar>();
   };

   double h;
  private:

   Matrix<_Scalar, Dynamic, 1> x0_;
   Matrix<_Scalar, Dynamic, 1> dx0_;
   Matrix<_Scalar, Dynamic, 1> val1, val2;
   
   Matrix<_Scalar, Dynamic, 1> u_;
   Matrix<_Scalar, Dynamic, 1> x_;

   Matrix<_Scalar, Dynamic, Dynamic> jac; // jacobian
   double eps;
   int N_x, 
       N_z;

   // update jacobians 
   int updateJ(       Matrix<_Scalar, Dynamic, Dynamic> &J,
                const Matrix<_Scalar, Dynamic, 1> y0, // value to differentiate
                const Matrix<_Scalar, Dynamic, 1> u0, // constant input
                const bool diffx 
              )
   { if(diffx)
       jac.resize(x.rows(), x.rows());
     else
       jac.resize(x.rows(),u0.rows());
     int nfev=0;
     const typename int n = x.size()*int(diffx) +u0.size()*int(!diffx);
     eps = DBL_EPSILON;//std::sqrt(((std::max)(0.0,DBL_EPSILON)));
     x_ = x;
     u_ = u0;
     Matrix<_Scalar, Dynamic, 1> delta;
     if(diffx)
     { val1.resize(x_.rows());
       val2.resize(x_.rows());
       delta.resize(x_.rows());
       delta=eps *x_.cwiseAbs();
     }
     else
     { delta.resize(u_.rows());
       delta=eps *u_.cwiseAbs();
     }
  
     for (int j = 0; j < n; ++j) 
     { if (delta(j) == 0.) 
         delta(j) = eps;
       if(diffx)
       { x_(j) += delta(j);
       }
       else
         u_(j) += delta(j);
       x_.head(N_x);
       x_.segment(N_x,N_z);
       _Functor::operator()(x_.head(N_x), 
                            u_, 
                            x_.segment(N_x,N_z), 
                            val1); nfev++;
       if(diffx)
         x_(j) -= 2*delta(j);
       else
         u_(j) -= 2*delta(j);
       _Functor::operator()(x_.head(N_x), 
                            u_, 
                            x_.segment(N_x, N_z), 
                            val2); nfev++;
       if(diffx)
         x_ = x;
       else
         u_ = u0;
       //std::cout<<"val1: "<< val1.transpose()<<std::endl;
       //std::cout<<"val2: "<< val2.transpose()<<std::endl;
       jac.col(j) = (val2-val1)/(2*delta(j));
     }
     J = jac;
     return nfev;
   };
};

// central differencing
//H(\omega)=2i\sum_{k=1}^{(N-1)/2}{c_k\sin(k\omega)}
//Magnitude responses for N = 3,5,7,9 are drawn below:
// third order point stencil
// A= [dF/dx1 dF/dx2 dF/dxn]
// B= [dF/du1 dF/du2 dF/dum]

#endif