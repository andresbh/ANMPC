#include "EKF.h"

template <class _Functor, typename _Scalar >
EKF<_Functor, _Scalar>::EKF(
     //const _Functor& F0,
     const Matrix<_Scalar, Dynamic, 1> x0,
     const Matrix<_Scalar, Dynamic, 1> y0,
     const Matrix<_Scalar, Dynamic, 1> u0,
     const double h0,
     const Matrix<double, Dynamic, Dynamic> C0,
     const Matrix<double, Dynamic, Dynamic> P0,
     const Matrix<double, Dynamic, Dynamic> Q0,
     const Matrix<double, Dynamic, Dynamic> R0
     )
 : EKF::LKF( x0, y0, P0, Q0, R0, MatrixXd::Zero(C0.cols(),C0.cols()), //A
                                 MatrixXd::Zero(C0.cols(),u0.rows()),     //B
                                 C0,                                      //C
                                 MatrixXd::Zero(C0.rows(),u0.rows())      //D
           )
  {   x = x0;
      h = h0;
      dx0_.resize(x.size());
      x0_.resize(x.size());
//      model = F0;
  }

template <class _Functor, typename _Scalar >
int EKF<_Functor,_Scalar>::operator()(
     const Matrix<_Scalar, Dynamic, 1> x0
     )
   {  x = x0;
      return 0;
   }

template <class _Functor, typename _Scalar >
void EKF<_Functor,_Scalar>::update(
     const Matrix<_Scalar, Dynamic, 1> y0, 
     const Matrix<_Scalar, Dynamic, 1> u0
    )
   {  //Matrix<double, Dynamic, 1> y00_=y0.cast<double>();
      //Matrix<double, Dynamic, 1> u00_=u0.cast<double>();

      updateJ(A, y0, u0, true); // continuous state space
      A = MatrixXd::Identity(A.rows(), A.cols()).cast<AD<double>>() +h*A;
      updateJ(B, y0, u0, false); // continuous state space
      B = h*B;
      x0_ = x;
      LKF::update(y0, u0);

//      Matrix<double, Dynamic, 1> x00_=x0_.cast<double>();
//      Matrix<double, Dynamic, 1> dx00_=dx0_.cast<double>();

      _Functor::operator()(x_, u_, dx0_);

      x = (x_ +dx0_*h).cast<_Scalar>();
      return;
   }

template <class _Functor, typename _Scalar >
int EKF<_Functor,_Scalar>::updateJ( // update jacobians 
           Matrix<_Scalar, Dynamic, Dynamic> &J,
     const Matrix<_Scalar, Dynamic, 1> y0, // value do differentiate
     const Matrix<_Scalar, Dynamic, 1> u0, // constant input
     const bool diffx 
    )
   { if(diffx)
       jac.resize(x.rows(), x.rows());
     else
       jac.resize(x.rows(),u0.rows());
     //using std::abs;
     // Local variables 
     // double delta;
     int nfev=0;
     const typename int n = x.size()*int(diffx) +u0.size()*int(!diffx);
     eps = std::sqrt(((std::max)(0.0,DBL_EPSILON)));
     x_ = x;
     u_ = u0;
     val1.resize(x_.rows());
     val2.resize(x_.rows());
     Matrix<_Scalar, Dynamic, 1> delta =eps *x_.cwiseAbs();
     
     for (int j = 0; j < n; ++j) 
     { if (delta(j) == 0.) 
         delta(j) = eps;
       if(diffx)
         x_(j) += delta(j);
       else
         u_(j) += delta(j);
       _Functor::operator()(x_, u_, val2); nfev++;
       if(diffx)
         x_(j) -= 2*delta(j);
       else
         u_(j) -= 2*delta(j);
       _Functor::operator()(x_, u_, val1); nfev++;
       if(diffx)
         x_ = x;
       else
         u_ = u0;
       jac.col(j) = (val2-val1)/(2*delta(j));
     }
     J = jac;
     return nfev;
   }

//The explicit instantiation part
//template class EKF<model, double>;
template class EKF<model, AD<double>>;
template class LKF< COSTF::ADdouble >;