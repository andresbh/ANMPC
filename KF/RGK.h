#pragma once
#ifndef RGK_H
#define RGK_H
#define EIGEN_FFTW_DEFAULT
#include <Eigen/Eigen>
using namespace Eigen;

class RGK
{
public:
  RGK(Matrix<double, Dynamic, 1> x0);
  ~RGK(void);

  protected:
   Matrix<double, Dynamic, 1> f(Matrix<double, Dynamic, Dynamic> x, Matrix<double, Dynamic, Dynamic> y);

  public:
};

#endif
