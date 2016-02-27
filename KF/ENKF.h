#pragma once
#include "LKF.h"

template <class _Functor, typename _Scalar>
class ENKF : public LKF<_Functor, _Scalar>
{
public:
  ENKF(void);
  ~ENKF(void);
};
