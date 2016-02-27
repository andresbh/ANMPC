# include "NMPCFormulation.hpp"
NMPC::NMPC()
{ }

NMPC::NMPC(const int N_x0,       // number of states
           const int N_y0,       // number of outputs
           const int N_z0,       // number of controlled vars
           const int N_u0,       // number of inputs
           const int H_x0,       // prediction horizon
           const VectorXi ublk0, // input prediction blocking variable
           const double h0       // cycle time
           ) 
  : FG_info(N_x0, N_y0, N_z0, N_u0, H_x0, ublk0, h0)
{ int i, j, k;
  this->init_b = true;

  this->ok = true;
  this->h = h0;

  n = N_u0 *ublk0.size();
  m = N_z0;
  u_applied.resize(n, 1);
  u_applied.fill(0.0);
  u_l.resize(n);
  u_h.resize(n);

  for(i=0, k=0;i<N_u0; i++)
  { for(j=0; j<ublk0.size(); j++, k++)
    { u_l(k)=(0.0);
      u_h(k)=(1.0E1);
    }
  }
  // lower and upper limits for g
  g_l.resize(m);
  g_h.resize(m);
  for(i=0;i<m;i++)
  { g_l(i)=(0.0);
    g_h(i)=(2.0E1);
  }

  // retape operation sequence for each new u
  options += "Retape  false\n";            // retape
  //options += "Sparse  true  forward\n"; // use sparse jacobian and hesian with forward
  // turn off any printing
  options += "Integer print_level      0\n"; 
  options += "String  sb               yes\n";
  // maximum number of iterations
  options += "Integer max_iter         10\n";
  // maximum time of iterations in seconds
  options += "Numeric max_cpu_time     1\n";
  // approximate accuracy in first order necessary conditions;
  // see Mathematical Programming, Volume 106, Number 1, 
  // Pages 25-57, Equation (6)
  options += "Numeric tol            1e-7\n";
  // derivative testing
  //options += "String  derivative_test            second-order\n";
  options += "String  derivative_test            first-order\n";
  // maximum amount of random pertubation; e.g., 
  // when evaluation finite diff
  //options += "Numeric point_perturbation_radius  0e0\n";
  //options += "String hessian_approximation  limited-memory\n"; // inacurate solving wiht EKF
  //options += "String mu_strategy adaptive\n";

  solution = new CppAD::ipopt::solve_result<Dvector>();

  Matrix<ADdouble, Dynamic, 1> uapplied0_ = this->u_applied.cast<ADdouble>().head(N_u);
  
  this->EKFS->update(this->measurements, uapplied0_);
  
  CppAD::ipopt::solve<VectorXd, FG_info>(options, u_applied, u_l, u_h, g_l, g_h, *this, *solution);

  // Check some of the solution values
  this->ok = true;//solution->status == CppAD::ipopt::solve_result<Dvector>::success;
  //this->u_applied = solution->x;
}

NMPC::~NMPC()
{  
}

bool NMPC::update()
{ // place to return solution
  // std::cout<<"optimSolution: "<<solution->x<<std::endl;
  bool debug_b = false;
  if(!this->init_b)
  { u_applied = this->solution->x; // should read from last applied system input
    //std::cout<<u_applied<<std::endl<<this->solution->x<<std::endl;
  }
  else
    this->init_b = false;
  // update kalman filter
  Matrix<ADdouble, Dynamic, 1> uapplied0_ = this->u_applied.cast<ADdouble>().head(N_u);
  this->EKFS->update(this->measurements, uapplied0_);

  //std::cout<<"uapp: "<<u_applied.transpose()<<std::endl;
  CppAD::ipopt::solve<VectorXd, FG_info>(options, u_applied, u_l, u_h, g_l, g_h, *this, *solution);
  switch(solution->status)
  { case CppAD::ipopt::solve_result<Dvector>::success:
    case CppAD::ipopt::solve_result<Dvector>::maxiter_exceeded:
	  if (debug_b)
		std::cout<<"maxiter_exceeded"<<std::endl;
    case CppAD::ipopt::solve_result<Dvector>::stop_at_acceptable_point:
	  if (debug_b)
		std::cout<<"stop_at_acceptable_point"<<std::endl;
    case CppAD::ipopt::solve_result<Dvector>::unknown:
	  if (debug_b)
		std::cout<<"unknown"<<std::endl;
      u_applied = solution->x; 
      this->ok = true;
      break;
    default:
      this->ok = true;
      break;
  }

  return this->ok;
}

bool NMPC::solve()
{   // Check some of the solution values
  return this->ok;
}

// END C++
