#pragma once
#ifndef UKF_H
#define UKF_H

template <class _Functor>
class UKF: public _Functor 
{
private:

public:
  _Functor Functor;
  double h, lambda, c, csqrt;
  double alpha, beta;
  double ki;
  Matrix<double, Dynamic, 1> Wm, Wc;
  Matrix<double, Dynamic, Dynamic> X, Y, K, C, P, Q, R, A;
  Matrix<double, Dynamic, 1> x, dx0_;

  UKF(
    const Matrix<double, Dynamic, 1> x0,
    const Matrix<double, Dynamic, 1> y0,
    const Matrix<double, Dynamic, 1> u0,
    const double h0,
    const Matrix<double, Dynamic, Dynamic> C0,
    const Matrix<double, Dynamic, Dynamic> P0,
    const Matrix<double, Dynamic, Dynamic> Q0,
    const Matrix<double, Dynamic, Dynamic> R0
    ):x(x0),dx0_(x0),h(h0),alpha(1.0e-3),ki(0.0),beta(2.0),
      C(C0),P(P0),Q(Q0),R(R0)
  { int L = x0.size();
    
    this->X.resize(x.rows(),1+2*x.rows());
    this->X.fill(0.0);
    this->lambda = alpha*alpha*(L+ki)-L;
    this->c=x.size()+lambda;
    this->Wm.resize((2*x.size()+1),1);
    this->Wm(0)=lambda/c;
    this->Wm.tail(2*x.size()).fill(0.5/c);
    this->Wc=Wm;
    this->Wc(0)=Wc(0)+(1-alpha*alpha+beta);
    this->csqrt= sqrt(c);
    this->Y.resize(x.size(),1+2*x.size());
    cout<<"X: "<<X<<endl;
    cout<<"lambda: "<<lambda<<endl;
    cout<<"c: "<<c<<endl;
    cout<<"Wm: "<<Wm<<endl;
    cout<<"Wc: "<<Wc<<endl;
    cout<<"csqrt: "<<csqrt<<endl;
    cout<<"Y: "<<this->Y<<endl;
  };

  void update(
    const Matrix<double, Dynamic, 1> y0, 
    const Matrix<double, Dynamic, 1> u0
   )
  {  Matrix<double, Dynamic, 1> x0, x1, z1;
     Matrix<double, Dynamic, Dynamic> P1, P2, P12, X1, X2, Z1, Z2;
     sigmas(X,x);

     this->unscentedTrans(x1, X1, P1, X2,
       u0, X , Wm, Wc, x.size() , Q, &_Functor::processDynamics, Functor);
     this->unscentedTrans(z1, Z1, P2, Z2,
       u0, X1, Wm, Wc, y0.size(), R, &_Functor::processOutput, Functor);

     P12 = X2*Wc.diagonal()*Z2.transpose();
     K=P2.lu().solve(P12);
     x=x1+K*(y0-z1);
     P=P1-K*P12.transpose();
  };

  void
   sigmas(       Matrix<double, Dynamic, Dynamic> &X0,
           const Matrix<double, Dynamic, 1>        x0)
  { //cout<<"P: " <<P<<endl;
    this->A = this->P.llt().matrixU().transpose();
    this->A*=this->c;
    Matrix<double, Dynamic, Dynamic> Y_m;
    Y_m.setZero(A.cols(),A.rows());
    for(int i=0;i<x0.size();i++)
    { Y_m.col(i) = x0;
      //cout<<"x:"<<x<<endl;
    }
    cout<<"Y_m: "<<Y_m<<endl;
    cout<<"A: "<<A<<endl;

    X0<<x0, Y_m+this->A, Y_m-this->A;
  };

  void
   unscentedTrans(         Matrix<double, Dynamic, 1>       &y0,
                           Matrix<double, Dynamic, Dynamic> &Y0,
                           Matrix<double, Dynamic, Dynamic> &P0,
                           Matrix<double, Dynamic, Dynamic> &Y1,
                     const Matrix<double, Dynamic, 1>       u0,
                     const Matrix<double, Dynamic, Dynamic> X0,
                     const Matrix<double, Dynamic, Dynamic> Wm0,
                     const Matrix<double, Dynamic, Dynamic> Wc0,
                     const int L,
                     const Matrix<double, Dynamic, Dynamic> Q0, 
                     void (_Functor::*processFunction)(
                                            const Matrix<double, Dynamic, 1>, 
                                            const Matrix<double, Dynamic, 1>, 
                                            const Matrix<double, Dynamic, 1>,
                                                  Matrix<double, Dynamic, 1>&
                                           ),
                     _Functor& Functor
                     )
  { Matrix<double, Dynamic, 1> x0;
    Matrix<double, Dynamic, Dynamic> x0_m;
    x0.resize(x.rows(),1);
    x0.fill(0.0);
    x0_m.resize(x.rows(),1);
    x0_m.fill(0.0);
    cout<<"Y.rows() : "<<this->Y.rows()<<endl;
    cout<<"Y.cols() : "<<this->Y.cols()<<endl;
    cout<<"x.rows() : "<<this->x.rows()<<endl;
    for (int i=0;i<1+2*x.rows();i++)
    { (Functor.*processFunction)(x0, u0, x, dx0_);
      Y.col(i) = x0;
      cout<<"Wm(i): "<<Wm(i)<<endl;
      cout<<"Y.col(i): "<<Y.col(i)<<endl;
      cout<<"x0_m: "<<x0_m<<endl;
      cout<<"x0_m+Wm(i)*Y.col(i): "<<x0_m+Wm(i)*Y.col(i)<<endl;

      x0_m+=Wm(i)*Y.col(i);
    }

    y0.resize(x0_m.row(0).size());
    cout<<"x0_m.row(0) : "<<x0_m.row(0)<<endl;
    cout<<"y0 : "<<y0<<endl;
    y0<<x0_m.row(0);
    Y1.resize(Y.rows(),X0.cols());
    for(int i=0;i<Y.cols();i++)
    { Y1.col(i) = Y.col(i)-x0_m;
    }
    cout<<"(Y1* this->Wc.asDiagonal() *Y1.transpose()) : "<<(Y1* this->Wc.asDiagonal() *Y1.transpose())<<endl;
    cout<<"Q0 : "<<Q0<<endl;
    P0 = (Y1* this->Wc.asDiagonal() *Y1.transpose()) +Q0;
  };

};


#endif