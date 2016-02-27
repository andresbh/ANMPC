#pragma once
#ifndef __NMPCFormulation_HPP__
#define __NMPCFormulation_HPP__

#include <iostream>
#include "realModel.h"
#include "optFormulation.hpp"

class NMPC : public FG_info
{ 
  private:

  public:
    NMPC();
    // derived class part of constructor
    NMPC(const int N_x0, const int N_y0, const int N_z0,
         const int N_u0, const int H_x0, 
         const VectorXi ublk0, const double h0);
    // default destructor //
    ~NMPC();
    // update optimization formulation
    bool update();
    //solve optimization
    bool solve();
    CppAD::ipopt::solve_result<COSTF::Dvector> *solution;

    COSTF::Dvector  u_l; // opt value llim
    COSTF::Dvector  u_h; // opt value hlim

    COSTF::Dvector  g_l; // opt nlin low limit
    COSTF::Dvector  g_h; // opt nlin high limit

    // options 
    std::string options; // optimizer options
  
    int n, // size of optimization variable
        m; // size of constraint functions
    bool init_b;
    double h;

    bool ok;

};


#endif
