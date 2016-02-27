#pragma once
#ifndef __OPTFormulation_HPP__
#define __OPTFormulation_HPP__

#include "Cost.hpp"

class FG_info : public COSTF
{ private:
  public:
   FG_info();
   // derived class part of constructor
   FG_info(const int N_x0, const int N_y0, const int N_z0, const int N_u0, const int H_x0, const VectorXi ublk0, const double h0);
   void operator()(ADvector& fg, const ADvector& u);
};

#endif
