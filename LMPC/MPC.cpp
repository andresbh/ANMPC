//#include <stdlib.h>

#include <qpOASES.hpp>

#include <iostream>

#include <Eigen/Eigen>
#include "LTI.h"

using namespace Eigen;
using namespace LTI_Object;
USING_NAMESPACE_QPOASES


template <typename _Scalar>
class LMPC: public LTI<_Scalar>
{  
public:
  Options mpcOptions;
  SQProblem *mpcQP;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Grun;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Omega;
  Eigen::Matrix<_Scalar, Dynamic, 1> omega;

  Eigen::Matrix<_Scalar, Dynamic, 1> r_traj;
  
  LMPC(Eigen::Matrix<_Scalar, Dynamic, 1> y0,
       Eigen::Matrix<_Scalar, 1, Dynamic> A0,
       Eigen::Matrix<_Scalar, 1, Dynamic> B0
    ) : LMPC::LTI(y0, A0, B0)
  {  };
  
  LMPC(const Eigen::Matrix<_Scalar, Dynamic, 1>       y0,
       const Eigen::Matrix<_Scalar, Dynamic, Dynamic> A0,
       const Eigen::Matrix<_Scalar, Dynamic, Dynamic> B0,
       const Eigen::Matrix<_Scalar, Dynamic, Dynamic> C0,
       const Eigen::Matrix<_Scalar, Dynamic, Dynamic> D0
    ) : LMPC::LTI(y0, A0, B0, C0, D0)
  {
  };

  void init(
    const Eigen::Matrix<_Scalar, Dynamic, 1> x_est,
    const Eigen::Matrix<_Scalar, Dynamic, 1> u_last
    )
  {
    r_traj.resize(this->Hp*this->B.cols(), 1);
    r_traj.fill(0);

    Omega = this->F_con;
    this->append(Omega, this->W_con);
    
    omega = -F_con.block(0, 0, F_con.rows(), this->B.cols())*u_last + this->f_con;
    this->append(omega, this->w_con);

    r_traj -= this->Phi*x_est - this->Gamma*u_last;
    Grun = 2 * this->Theta.transpose()*this->Qq*r_traj;

    std::cout << "Omega: " << std::endl << Omega << std::endl;
    std::cout << "omega: " << std::endl << omega << std::endl;
    
    mpcOptions.setToMPC();
    mpcOptions.printLevel = PL_LOW;
    int na_l=this->H.rows();
    int nb_l = this->Omega.rows();
    
    mpcQP = new SQProblem(na_l, nb_l);

    mpcQP->setOptions(mpcOptions);
    //SolutionAnalysis myAnalysis;
    Eigen::Matrix<real_t, Eigen::Dynamic, 1> lbAa, lbb, ubb;
    lbb.resize(this->H.cols(), 1);
    lbb.fill(-1e32);
    ubb.resizeLike(lbb);
    ubb.fill(1e32);
    lbAa.resize(this->Omega.cols(),1);
    lbAa.fill(-1e32);
    int nWSR = 100;
    real_t cputime;
    std::cout << "H:" << std::endl << this->H << std::endl;
    std::cout << "Grun:" << std::endl << this->Grun << std::endl;
    std::cout << "Omega:" << std::endl << this->Omega << std::endl;
    std::cout << "lbb:" << std::endl << lbb << std::endl;
    std::cout << "ubb:" << std::endl << ubb << std::endl;
    std::cout << "lbAa:" << std::endl << lbAa << std::endl;
    std::cout << "omega:" << std::endl << omega << std::endl;
    mpcQP->init(this->H.data(), this->Grun.data(), this->Omega.data(),
      lbb.data(), ubb.data(),
      lbAa.data(), this->omega.data(), nWSR, &cputime);
  };
  
  ~LMPC() 
  {
    delete mpcQP;
  };

  // qp solution --------------------------------------------
  double qpSolve
    (
      const Eigen::Matrix<_Scalar, Dynamic, Dynamic> x_est, // state estimate
      const Eigen::Matrix<_Scalar, Dynamic, Dynamic> u_last, // last input
      const Eigen::Matrix<_Scalar, Dynamic, Dynamic> du_old, // last input change
      const Eigen::Matrix<_Scalar, Dynamic, Dynamic> r,  // reference
            Eigen::Matrix<_Scalar, Dynamic, Dynamic> &u  // output
      )
  {
    for (int ii=0;ii < r.rows();ii++)
    {
      r_traj.block(ii, 0, ii*this->Hp, 1).fill(r(ii));
    }

    r_traj -= this->Phi*x_est - this->Gamma*u_last;
    Grun = 2 * this->Theta.transpose()*this->Qq*r_traj;

    Omega = this->F_con;
    this->append(Omega, this->W_con);
    
    omega = -F_con.block(0, 0, F_con.rows(), this->B.cols())*u_last + this->f_con;
    this->append(omega, this->w_con);
    std::cout << "Omega: " << std::endl << Omega << std::endl;
    std::cout << "omega: " << std::endl << omega << std::endl;
    
    Eigen::Matrix<real_t, Eigen::Dynamic, 1> lbAa, lbb, ubb;
    lbb.resize(this->H.cols(), 1);
    lbb.fill(-1e32);
    ubb.resizeLike(lbb);
    ubb.fill(1e32);
    lbAa.resize(this->Omega.cols(), 1);
    lbAa.fill(-1e32);

    int nWSR = 10;
    double cputime;
    mpcQP->hotstart(this->H.data(), this->Grun.data(), this->Omega.data(),
      lbb.data(), ubb.data(),
      lbAa.data(), this->omega.data(), nWSR, 0);
    std::cout << "u: " << u << std::endl;
    mpcQP->getPrimalSolution(u.data());
    std::cout << "u: " << u << std::endl;
    return mpcQP->getObjVal(); // value of cost function
  };
};
 
int main( )
{ //USING_NAMESPACE_QPOASES

  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> Hh, Aa, gg, lbAa, ubAa, lbb, ubb;
  
  //Eigen::Matrix<real_t, 1, Eigen::Dynamic> A0(1, 1); A0 << -0.32870;
  //Eigen::Matrix<real_t, 1, Eigen::Dynamic> B0(1, 1); B0 << 0.6065;
  
  LMPC<real_t> *LMPC_O;
  Eigen::VectorXd x,y;
  y.resize(2); y.setZero();
  x.resize(1); x.setZero();
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> A0(4, 4); 
  A0<< 0.94301  , 0.00000 , 0.07564 ,  0.00000,
       0.00000  , 0.97323 , 0.00000 ,  0.03857,
       0.00000  , 0.00000 , 0.92209 ,  0.00000,
       0.00000  , 0.00000 , 0.00000 ,  0.96090;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> B0(4, 2);
  B0<< 0.08663 ,  0.00903,
       0.00457 ,  0.10844,
       0.00000 ,  0.22409,
       0.22953 ,  0.00000;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> C0(2, 4); 
  C0<< 0.50000 ,  0.00000 , 0.00000 , 0.00000,
       0.00000 ,  0.50000 , 0.00000 , 0.00000;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> D0(2, 2); 
  D0<< 0 ,  0,
       0 ,  0;

  Eigen::Matrix<double, Dynamic, Dynamic> u0, x0, Du0, r;
  u0.resize(2, 3);
  u0 << 0, 0, 0, 0, 0, 0;
  x0.resize(A0.cols(), 1);
  x0 << 1, 0, 0, 0;
  r.resize(C0.rows(), 1);
  r << 1, 0;

  LMPC_O = new LMPC<real_t>( y, A0, B0, C0, D0);
  
  Du0.resize(LMPC_O->Theta.cols(), 1);
  Du0.setZero();

  LMPC_O->createPredictionMatrices(5, 3);

  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> Q(2, 2);
  Q << 4, 0,
       0, 1;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> R(2, 2);
  R << 1, 0,
       0, 1;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> umax(2, 1), umin(2, 1), dumax(2, 1) , dumin(2, 1);
  umax << 1, 2;umin << -3, -4;dumax << 5, 6;dumin << -7, -8;
  LMPC_O->createOptimizationMatrices(Q,R,umax,umin,dumax,dumin);

// not working LMPC_O->prediction(u0, Du0, x0);
  LMPC_O->init(x0, u0.block(0, 0, u0.rows(), 1));
  LMPC_O->qpSolve(x0, u0.col(0), Du0, r, u0);

	// Setup data of first QP...
  int na_l = 2, nb_l = 1;
  Hh.resize(na_l, na_l);
  Hh << 1.0, 0.0, 0.0, 0.5;
  Aa.resize( nb_l, na_l);
  Aa << 1.0, 1.0;
  gg.resize(na_l,1);
  gg << 1.5, 1.0;
  lbb.resize(na_l, 1);
  lbb << 0.5, -2.0;
  ubb.resize(na_l, 1);
  ubb << 5.0, 2.0;
  lbAa.resize(nb_l,1);
  lbAa << -1.0;
  ubAa.resize(nb_l, 1);
  ubAa << 2.0;

  Options myOptions;
  myOptions.setToMPC();
  myOptions.printLevel = PL_LOW;

  SQProblem myExample(na_l,nb_l);

  myExample.setOptions(myOptions);
  SolutionAnalysis myAnalysis;

  int nWSR = 10;
  real_t cputime;
  myExample.init(Hh.data(), gg.data(), Aa.data(), 
    lbb.data(), ubb.data(),
    lbAa.data(), ubAa.data(), nWSR, &cputime);
  
  // Get and print solution of first QP.
  real_t xOpt[2];
  real_t yOpt[2 + 1];
  myExample.getPrimalSolution(xOpt);
  myExample.getDualSolution(yOpt);
  printf("\nxOpt = [ %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n",
    xOpt[0], xOpt[1], yOpt[0], yOpt[1], yOpt[2], myExample.getObjVal());

  real_t H[2 * 2] = { 1.0, 0.0, 0.0, 0.5 };
  real_t A[1 * 2] = { 1.0, 1.0 };
  real_t g[2] = { 1.5, 1.0 };
  real_t lb[2] = { 0.5, -2.0 };
  real_t ub[2] = { 5.0, 2.0 };
  real_t lbA[1] = { -1.0 };
  real_t ubA[1] = { 2.0 };

  // Setting up SQProblem object and solution analyser.
  SQProblem example(2, 1);
  example.setOptions(myOptions);
  example.init(H, g, A, lb, ub, lbA, ubA, nWSR, &cputime);
  
  /* Get and print solution of first QP. */
  example.getPrimalSolution(xOpt);
  example.getDualSolution(yOpt);
  printf("\nxOpt = [ %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n",
    xOpt[0], xOpt[1], yOpt[0], yOpt[1], yOpt[2], example.getObjVal());

  real_t maxKktViolation = myAnalysis.getKktViolation(&myExample);
  printf("maxKktViolation: %e\n", maxKktViolation);
  
  maxKktViolation = myAnalysis.getKktViolation(&example);
  printf("maxKktViolation: %e\n", maxKktViolation);

  Hh << 1.0, 0.50, 0.50, 0.5;
  Aa << 1.0, 5.0;
  gg << 1.0, 1.5;
  lbb << 0.0, -1.0;
  ubb << 5.0, -0.50;
  lbAa << -2.0;
  ubAa << 1.0;
  
  nWSR = 10;
  myExample.hotstart(Hh.data(), gg.data(), Aa.data(), lbb.data(), ubb.data(),
    lbAa.data(), ubAa.data(), nWSR, 0);

  /* Setup data of second QP. */
  real_t H_new[2 * 2] = { 1.0, 0.5, 0.5, 0.5 };
  real_t A_new[1 * 2] = { 1.0, 5.0 };
  real_t g_new[2] = { 1.0, 1.5 };
  real_t lb_new[2] = { 0.0, -1.0 };
  real_t ub_new[2] = { 5.0, -0.5 };
  real_t lbA_new[1] = { -2.0 };
  real_t ubA_new[1] = { 1.0 };

  nWSR = 10;
  example.hotstart(H_new, g_new, A_new, lb_new, ub_new, lbA_new, ubA_new, nWSR, 0);

  myExample.getPrimalSolution(xOpt);
  myExample.getDualSolution(yOpt);
  printf("\nxOpt = [ %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n",
    xOpt[0], xOpt[1], yOpt[0], yOpt[1], yOpt[2], example.getObjVal());

  example.getPrimalSolution(xOpt);
  example.getDualSolution(yOpt);
  printf("\nxOpt = [ %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n",
    xOpt[0], xOpt[1], yOpt[0], yOpt[1], yOpt[2], example.getObjVal());


  maxKktViolation = myAnalysis.getKktViolation(&myExample);
  printf("maxKktViolation: %e\n", maxKktViolation);

  maxKktViolation = myAnalysis.getKktViolation(&example);
  printf("maxKktViolation: %e\n", maxKktViolation);

  //  ------------ VARIANCE-COVARIANCE EVALUATION --------------------

  real_t *Var = new real_t[5 * 5];
  real_t *Primal_Dual_Var = new real_t[5 * 5];
  real_t *myPrimal_Dual_Var = new real_t[5 * 5];

  int run1, run2;
  for (run1 = 0; run1 < 5 * 5; run1++)
    Var[run1] = 0.0;

  Var[0] = 1.0;
  Var[6] = 1.0;

  //                  (  1   0   0   0   0   )
  //                  (  0   1   0   0   0   )
  //     Var     =    (  0   0   0   0   0   )
  //                  (  0   0   0   0   0   )
  //                  (  0   0   0   0   0   )

  myAnalysis.getVarianceCovariance(&myExample, Var, myPrimal_Dual_Var);
  myAnalysis.getVarianceCovariance(&example, Var, Primal_Dual_Var);

  printf("\nPrimal_Dual_VAR = \n");
  for (run1 = 0; run1 < 5; run1++) {
    for (run2 = 0; run2 < 5; run2++) {
      printf(" %10f", myPrimal_Dual_Var[run1 * 5 + run2]);
    }
    printf("\n");
  }

  printf("\nPrimal_Dual_VAR = \n");
  for (run1 = 0; run1 < 5; run1++) {
    for (run2 = 0; run2 < 5; run2++) {
      printf(" %10f", Primal_Dual_Var[run1 * 5 + run2]);
    }
    printf("\n");
  }

  delete[] Primal_Dual_Var;
  delete[] Var;

  return 0;
}
