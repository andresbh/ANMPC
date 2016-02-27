# include "optFormulation.hpp"

FG_info::FG_info()
{ }

FG_info::FG_info( const int N_x0, 
                  const int N_y0, 
                  const int N_z0, 
                  const int N_u0, 
                  const int H_x0, 
                  const VectorXi ublk0, 
                  const double h0)
:  FG_info::COSTF(N_x0, N_y0, N_z0, N_u0, H_x0, ublk0.size(), ublk0, h0)
{ }

// Evaluation of the objective f(x), and constraints g(x)
// using an Algorithmic Differentiation (AD) class.
void FG_info::operator()(ADvector& fg, const ADvector& u)
{ int i=0, k=0;
  // cost function
  AD<double> val = costfct( u);
  fg(0) = val;

  // g_1 (x)
  fg.segment(1,N_z) = measurements.head(N_z);

  return;
}