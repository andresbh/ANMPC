#pragma once

#include <qpOASES.hpp>
#include <iostream>

#include <Eigen/Eigen>
#include "LTI.h"

using namespace Eigen;
using namespace LTI_Object;
USING_NAMESPACE_QPOASES

template <typename _Scalar>
class LMPC : public LTI<_Scalar>
{
public:
  typedef typename Eigen::Matrix<_Scalar, Dynamic, Dynamic> MatrixS; 
  Options mpcOptions;
  SQProblem *mpcQP;
  int na, ma, nb, mb, nc, mc;
  int py; // Number of measured outputs for system
  int pz; // Number of controlled outputs for the system
  int pc; // Number of constrained outputs for the system

  int pze, ne, pye, me;

  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Grun;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Omega;
  Eigen::Matrix<_Scalar, Dynamic, 1> omega;
  Eigen::Matrix<_Scalar, Dynamic, 1> r_traj;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Ad;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Bd;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Cyd;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Czd;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Dzd;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Ccd;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Dcd;
  
  Eigen::Matrix<int, Dynamic, 1> ublk;
  Eigen::Matrix<int, Dynamic, 1> duSampInclude;
  Eigen::Matrix<int, Dynamic, 1> zblk;
  Eigen::Matrix<int, Dynamic, 1> zSampInclude;

  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Ade;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Bde;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic>Cyde;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic>Czde;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic>Dzde;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic>Ccde;
  Eigen::Matrix<_Scalar, Dynamic, Dynamic>Dcde;

  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Psi; // prediction matrix
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Gamma; // prediction matrix
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Theta; // prediction matrix
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Psi_c; // constraint prediction matrix
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Gamma_c; // constraint prediction matrix
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Theta_c; // constraint prediction matrix

  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Qq; // optimization matrix
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Rr; // optimization matrix
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> H; // optimization matrix
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> G; // optimization matrix
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Tc; // predi matrix
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> Gam_con; // optimization constraints
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> gam_con; // optimization constraints
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> F_con; // optimization constraints
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> f_con; // optimization constraints
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> W_con; // optimization constraints
  Eigen::Matrix<_Scalar, Dynamic, Dynamic> w_con; // optimization constraints

  int Hp,Hu,Hw;

  // constructor -----------------------
  LMPC(Eigen::Matrix<_Scalar, Dynamic, 1> y0,
    Eigen::Matrix<_Scalar, 1, Dynamic> A0,
    Eigen::Matrix<_Scalar, 1, Dynamic> B0
    ) : LMPC::LTI(y0, A0, B0)
  {  };
  
  // constructor -----------------------
  /*
    x(k + 1) = Ax(k) + Bu(k) // states
    y(k) = Cyx(k)            // measured output
    z(k) = Czx(k) + Dzu(k)   // controlled output
    zc(k) = Ccx(k) + Dcu(k)  // constrained output

    ∆umin ≤ ∆u(k) ≤ ∆umax
    umin ≤ u(k) ≤ umax
    zmin ≤ zc(k) ≤ zmax
  */

  LMPC(const Eigen::Matrix<_Scalar, Dynamic, 1>    y0,
    const Eigen::Matrix<_Scalar, Dynamic, Dynamic> Ad,
    const Eigen::Matrix<_Scalar, Dynamic, Dynamic> Bd,
    const Eigen::Matrix<_Scalar, Dynamic, Dynamic> Cyd,
    const Eigen::Matrix<_Scalar, Dynamic, Dynamic> Czd,
    const Eigen::Matrix<_Scalar, Dynamic, Dynamic> Dzd,
    const Eigen::Matrix<_Scalar, Dynamic, Dynamic> Ccd,
    const Eigen::Matrix<_Scalar, Dynamic, Dynamic> Dcd,
    const Eigen::Matrix<int, Dynamic, 1> zblk,
    const Eigen::Matrix<int, Dynamic, 1> ublk,
    int Hp, int Hu, int Hw
    ) : LMPC::LTI(y0, Ad, Bd, Cyd, Dzd),
        Ad(Ad),  Bd(Bd),  Cyd(Cyd),  Czd(Czd),  
        Dzd(Dzd),  Ccd(Ccd),  Dcd(Dcd),
        zblk(zblk), ublk(ublk),
        Hp(Hp), Hu(Hu), Hw(Hw)
  {
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> temp;
    na = Ad.rows();  ma = Ad.cols();
    nb = Bd.rows();  mb = Bd.cols();
    nc = Cyd.rows(); mc = Cyd.cols();
    
    py = nc;
    pc = Ccd.cols();
    pz = Czd.rows();

    Ade.resizeLike(A);// [Ad Bd; zeros(mb, na) eye(mb)];
    Ade = A;
    appendh(Ade, Bd);
    temp.resize(mb,na); temp.setZero();
    appendh(temp, MatrixS::Identity(mb, mb));
    appendv(Ade,temp);

    Bde.resize(nb, mb); //[Bd; zeros(mb, mb)];
    Bde = Bd;
    appendv(Bde, MatrixS::Zero(mb, mb));

    Cyde.resizeLike(Cyd); // = [Cyd zeros(py, mb)];
    Cyde = Cyd;
    appendh(Cyde, MatrixS::Zero(py, mb));

    Czde.resizeLike(Czd); // = [Czd zeros(pz, mb)];
    Czde = Czd;
    appendh(Czde, MatrixS::Zero(pz, mb));

    Dzde.resizeLike(Dzd);
    Dzde = Dzd;

    Ccde.resizeLike(Ccd); //[Ccd zeros(pc, mb)];
    Ccde = Ccd;
    appendh(Ccde, MatrixS::Zero(pc,mb));
    
    Dcde.resizeLike(Dcd);
    Dcde = Dcd;

    if (py > mb)
    {
      blkdiag(Ade, MatrixS::Identity(py - mb, py - mb));
      appendv(Bde, MatrixS::Zero(py - mb, mb));
      
      Cyde.resizeLike(Cyd); Cyde = Cyd; // [Cyd [zeros(mb, py);[zeros(py - mb, mb) eye(py - mb)]]];
      temp = appendvr(MatrixS::Zero(mb, py), appendhr(MatrixS::Zero(mb, py), MatrixS::Zero(py - mb, py - mb)));
      appendv(Cyde,temp); 

      appendh(Czde, MatrixS::Zero(pz, py - mb));
      appendh(Ccde, MatrixS::Zero(pc, py - mb));
    }

    pze = Czde.rows();
    pye = Cyde.rows();
    ne  = Bde.rows();
    me = Bde.cols();

  };

  // init -----------------------
  void init(
    const Eigen::Matrix<_Scalar, Dynamic, 1> x_est,
    const Eigen::Matrix<_Scalar, Dynamic, 1> u_last
    )
  {
    r_traj.resize(this->Hp*this->B.cols(), 1);
    r_traj.fill(0);

    Omega = this->F_con;
    this->appendv(Omega, this->W_con);

    omega = -F_con.block(0, 0, F_con.rows(), this->B.cols())*u_last + this->f_con;
    this->append<_Scalar>(omega, this->w_con);

    r_traj -= this->Psi*x_est - this->Gamma*u_last;
    Grun = 2 * this->Theta.transpose()*this->Qq*r_traj;

    std::cout << "Omega: " << std::endl << Omega << std::endl;
    std::cout << "omega: " << std::endl << omega << std::endl;

    mpcOptions.setToMPC();
    mpcOptions.printLevel = PL_LOW;
    int na_l = this->H.rows();
    int nb_l = this->Omega.rows();

    mpcQP = new SQProblem(na_l, nb_l);

    mpcQP->setOptions(mpcOptions);
    //SolutionAnalysis myAnalysis;
    Eigen::Matrix<real_t, Eigen::Dynamic, 1> lbAa, lbb, ubb;
    lbb.resize(this->H.cols(), 1);
    lbb.fill(-1e32);
    ubb.resizeLike(lbb);
    ubb.fill(1e32);
    lbAa.resize(this->Omega.cols(), 1);
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

  // destructor -----------------------
  ~LMPC()
  {
    delete mpcQP;
  };

  // createBlockingMatricesMatrices -----------------------
  void createBlockingMatrices() 
  {
    int ublklast = ublk.rows();
    int zblklast = zblk.rows();

    if (zblklast >= Hp - 1)
      zblklast = Hp - 1;

    if (ublklast >= Hu - 1)
      ublklast = Hu - 1;

    //zSample
    Eigen::Matrix<int, Dynamic, 1> temp;
    temp.resize(1, 1); temp(0, 0) = Hw + 1;
    append<int>(temp, zblk.head(zblklast));
    append<int>(temp, zblk(zblklast - 1)*MatrixXi::Ones(Hp - zblk.rows() - 1, 1));
    zSampInclude.resizeLike(temp);
    zSampInclude = cumsum<int>(temp);

    //duSample
    temp.resize(1, 1); temp(0, 0) = 0;
    append<int>(temp, ublk.head(ublklast));

    append<int>(temp, ublk(ublklast - 1)*MatrixXi::Ones(Hu - ublk.rows() - 1, 1));
    duSampInclude.resizeLike(temp);
    duSampInclude = cumsum<int>(temp);

  }
  
  // createPredictionMatrices -----------------------
  void createPredictionMatrices() // prediction horizon, control horizon
  {
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> T, Tc;
    Eigen::Matrix<int, Dynamic, Dynamic> zRowInclude;
    T.resizeLike(Czde); 
    Tc.resizeLike(Ccde);

    Eigen::MatrixPower<Eigen::Matrix<_Scalar, Dynamic, Dynamic>> Apow(this->Ade);
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> temp;
    temp.resize(ne, Gamma.cols());

    for (int kk = 1;kk < zSampInclude.tail(1)(0)+2; kk++)
    {
      if(kk==1)
      {
        T = Czde;
        Tc = Ccde;
        Psi.resizeLike(Ade); Psi = Ade;
        Gamma.resizeLike(Bde);Gamma = Bde;
        Theta.resizeLike(Bde);Theta = Bde;
        appendh(Theta, MatrixS::Zero(ne, me*this->duSampInclude.bottomRows(1)(0)));
      }
      else
      {
        this->blkdiag(T, this->Czde);
        this->blkdiag(Tc, this->Ccde);
        //populate Gamma
        appendv(Gamma, (Apow(kk - 1)*Bde + Gamma.block(ne*(kk - 2), 0, ne, Gamma.cols())));
        
        //populate Theta
        temp = Gamma.block(ne*(kk - 1), 0, ne, Gamma.cols());
        appendh(temp, Theta.block(ne*(kk - 2), 0, ne, me*duSampInclude.tail(1)(0)));
        appendv(Theta, temp);

        //populate Psi
        this->appendv(Psi, Apow(kk));
        // line 515
      }
    }

    Theta_c = Tc*Theta;
    Gamma_c = Tc*Gamma;
    Psi_c = Tc*Psi;

    Theta = T*Theta;
    Gamma = T*Gamma;
    Psi = T*Psi;
    
    // 604
    zRowInclude.resize(pze, zSampInclude.rows());
    zRowInclude.setOnes();
    std::cout << zRowInclude << std::endl<< std::endl;
    std::cout << zSampInclude << std::endl;
    zRowInclude = pze*zRowInclude.transpose()*zSampInclude;
    std::cout << zRowInclude << std::endl;
    //  zSampInclude;
    //zRowInclude;
  };


  // createOptimizationMatrices -----------------------
  void createOptimizationMatrices(
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> Q,
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> R,
    Eigen::Matrix<_Scalar, Dynamic, 1> umax,
    Eigen::Matrix<_Scalar, Dynamic, 1> umin,
    Eigen::Matrix<_Scalar, Dynamic, 1> dumax,
    Eigen::Matrix<_Scalar, Dynamic, 1> dumin,
    Eigen::Matrix<_Scalar, Dynamic, 1> zmax,
    Eigen::Matrix<_Scalar, Dynamic, 1> zmin
    )
  {
    int nbr_constr_ublk, nbr_constr_dublk;
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> temp;
    Eigen::Matrix<_Scalar, Dynamic, Dynamic> uf_block, uw_block, bu_block, au_block, z_block, gz_block, qwe(2,1), qwef(2, 1), qwe2(2, 1), qwew(2, 1);
    uf_block.resize(2 * B.cols(), B.cols()); uf_block.setZero();
    uw_block.resize(2 * B.cols(), B.cols()); uw_block.setZero();
    bu_block.resize(2 * B.cols(), 1); bu_block.setZero();
    au_block.resize(2 * B.cols(), 1); au_block.setZero();

    Qq.resizeLike(Q);
    Qq = Q;
    Rr.resizeLike(R);
    Rr = R;
    
    for (int ii = 0; ii < me; ii++) // umin or umax == inf not possible
    {
      qwef(0, 0) = 1;qwef(1, 0) = -1;
      uf_block.block(Bde.cols() * ii, ii, Bde.cols(), 1) = qwef;
      qwe2(0, 0) = umax(ii); qwe2(1, 0) = -umin(ii);
      bu_block.block(Bde.cols() * ii, 0, Bde.cols(), 1) = qwe2;

      qwe2(0, 0) = dumax(ii); qwe2(1, 0) = -dumin(ii);
      au_block.block(Bde.cols() * ii, 0, Bde.cols(), 1) = qwe2;

      qwew(0, 0) = 1;qwew(1, 0) = -1;
      qwe2(0, 0) = dumax(ii); qwe2(1, 0) = -dumin(ii);
      uw_block.block(Bde.cols() * ii, ii, Bde.cols(), 1) = qwef;
    }

    for (int ii = 0; ii < pc; ii++) // zmin or zmax
    {
      qwe(0, 0) = 1.0;qwe(1, 0) = -1.0;
      qwe2(0, 0) = zmax(ii);qwe2(1, 0) = -zmin(ii);
      if (ii == 0)
      {
        z_block.resizeLike(qwe);
        z_block=qwe;

        gz_block.resizeLike(qwe2);
        gz_block = qwe2;
      }
      else
      {
        blkdiag(z_block, qwe);
        appendv(gz_block, qwe2);
      }
    }
    // line 530
    f_con.resize(this->Hu * 2 * this->B.cols(), 1);
    f_con.setZero();

    nbr_constr_ublk = uf_block.rows();
    nbr_constr_dublk = uw_block.rows();

    for (int ii = 1; ii < zSampInclude.tail(1)(0); ii++)
    {
      blkdiag(Qq,Q);
    }

    for (int ii = 1; ii < duSampInclude.tail(1)(0); ii++)
    {
      blkdiag(Rr, R);
    }
    
    F_con.resizeLike(uf_block);
    F_con = uf_block;

    f_con.resizeLike(bu_block);
    f_con = bu_block;

    Gam_con.resizeLike(z_block);
    Gam_con = z_block;

    gam_con.resizeLike(gz_block);
    gam_con = gz_block;

    W_con.resizeLike(uw_block);
    W_con = uw_block;
    appendh(W_con, MatrixS::Zero(nbr_constr_ublk, me*(duSampInclude.bottomRows(1)(0) + 1) - me));

    w_con.resizeLike(au_block);
    w_con = au_block;

    this->appendh(F_con, 
      MatrixS::Zero(nbr_constr_ublk, me*(duSampInclude.bottomRows(1)(0)+1) - me));

    for (int ii = 1; ii < this->duSampInclude.bottomRows(1)(0)+1; ii++)
    {
      blkdiag(Rr, R);

      if (ii < zSampInclude.tail(1)(0) + 3)
      { 
        temp.resize(nbr_constr_ublk, ii*me);
        temp.setZero();
        appendh(temp, uf_block);
        appendh(temp, MatrixS::Zero(nbr_constr_ublk,me*(duSampInclude.bottomRows(1)(0) + 1) - (ii + 1)*me));
        appendv(F_con, temp);

        appendv(f_con, bu_block);

        temp.resize(nbr_constr_dublk, ii*me);
        temp.setZero();
        appendh(temp, uw_block);
        appendh(temp, MatrixS::Zero(nbr_constr_dublk, me*(duSampInclude.bottomRows(1)(0) + 1) - (ii + 1)*me));
        appendv(W_con, temp);

        appendv(w_con, au_block);
      }
    }

    for (int ii = this->duSampInclude.bottomRows(1)(0); ii >0 ; --ii)
    { //548
      F_con.block(0, (ii - 1)*me,F_con.rows(),me) += F_con.block(0, (ii)*me, F_con.rows(), me);
    }
    std::cout << "Rr" << std::endl;
    std::cout << Rr.rows()<<","<< Rr.cols() << std::endl;
    std::cout << "Qq" << std::endl;
    std::cout << Qq.rows() << "," << Qq.cols() << std::endl;
    std::cout << "Theta" << std::endl;
    std::cout << Theta.rows() << "," << Theta.cols() << std::endl;
    // 661
    H.resize(Rr.rows(), Rr.cols());
    H = Theta.transpose()*Qq*Theta + Rr;
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
    Z0.resize(this->Gamma.rows(), 1);

    Z0 = this->Psi*x0 + this->Gamma*u0.col(0) + this->Theta*Du0;

    //std::cout << "Z0:" << std::endl << Z0 << std::endl;

    return Z0;
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
    for (int ii = 0;ii < r.rows();ii++)
    {
      r_traj.block(ii, 0, ii*this->Hp, 1).fill(r(ii));
    }

    r_traj -= this->Psi*x_est - this->Gamma*u_last;
    Grun = 2 * this->Theta.transpose()*this->Qq*r_traj;

    Omega = this->F_con;
    appendv(Omega, this->W_con);

    omega = -F_con.block(0, 0, F_con.rows(), this->B.cols())*u_last + this->f_con;
    this->append<_Scalar>(omega, this->w_con);
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