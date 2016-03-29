#include "LMPC.h"
 
int main( )
{ Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> Hh, Aa, gg, lbAa, ubAa, lbb, ubb;

  //Eigen::Matrix<real_t, 1, Eigen::Dynamic> A0(1, 1); A0 << -0.32870;
  //Eigen::Matrix<real_t, 1, Eigen::Dynamic> B0(1, 1); B0 << 0.6065;

  LMPC<real_t> *LMPC_O;
  Eigen::VectorXd x,y;
  y.resize(2); y.setZero();
  x.resize(1); x.setZero();
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> Ad(4, 4); 
  Ad<< 0.94301  , 0.00000 , 0.07564 ,  0.00000,
       0.00000  , 0.97323 , 0.00000 ,  0.03857,
       0.00000  , 0.00000 , 0.92209 ,  0.00000,
       0.00000  , 0.00000 , 0.00000 ,  0.96090;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> Bd(4, 2);
  Bd<< 0.08663 ,  0.00903,
       0.00457 ,  0.10844,
       0.00000 ,  0.22409,
       0.22953 ,  0.00000;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> Cyd(2, 4); 
  Cyd<< 0.50000 ,  0.00000 , 0.00000 , 0.00000,
        0.00000 ,  0.50000 , 0.00000 , 0.00000;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> Dzd(2, 2); 
  Dzd<< 0 ,  0,
        0 ,  0;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> Czd(2, 4);
  Czd = Cyd;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> Ccd(4, 4);
  Ccd.setIdentity() *= 0.5;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> Dcd(4, 2);
  Dcd.setZero();

  Eigen::Matrix<double, Dynamic, Dynamic> u0, x0, Du0, r;
  Eigen::Matrix<int, Dynamic,1>  zblk, ublk;

  u0.resize(2, 3);
  u0 << 0, 0, 0, 0, 0, 0;
  x0.resize(Ad.cols(), 1);
  x0 << 1, 0, 0, 0;
  r.resize(Czd.rows(), 1);
  r << 1, 0;
  zblk.resize(2, 1);
  zblk << 1,2;
  ublk.resize(1, 1);
  ublk << 1;
  int Hp = 5, Hu = 2, Hw = 0;
  LMPC_O = new LMPC<real_t>(y, Ad, Bd, Cyd, Czd, Dzd, Ccd, Dcd, zblk, ublk, Hp, Hu, Hw);

  Du0.resize(LMPC_O->Theta.cols(), 1);
  Du0.setZero();

  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> Q(2, 2);
  Q << 4, 0,
       0, 1;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> R(2, 2);
  R << 1, 0,
       0, 1;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> zmin(4, 1), zmax(4, 1), umax(2, 1), umin(2, 1), dumax(2, 1), dumin(2, 1);
  umax << 1, 2;umin << -3, -4;dumax << 5, 6;dumin << -7, -8;
  zmax << 1E10, 2E10, 3E10, 4E10;zmin << -5E10, -6E10, -7E10, -8E10;
  
  LMPC_O->createBlockingMatrices();
  LMPC_O->createPredictionMatrices();
  LMPC_O->createOptimizationMatrices(Q,R,umax,umin,dumax,dumin,zmax,zmin);
  LMPC_O->prediction(u0, Du0, x0);
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
