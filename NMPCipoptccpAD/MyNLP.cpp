// MyNLP.c

#include "MyNLP.hpp"
using namespace Ipopt;

/* Constructor. */
NLMPC::NLMPC(const int N_x0, const int N_y0, const int N_u0, const int H_x0, const int H_u0, const int *ublk0 ):
  NLMPC::COSTF(N_x0, N_y0, N_u0, H_x0, H_u0, ublk0)
{  u_applied.resize(N_u0,1);
}

NLMPC::~NLMPC()
{}

bool NLMPC::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,Index& nnz_h_lag, IndexStyleEnum& index_style)
{ n = this->N_u*this->H_u;     // to input q flow at different sample times q[0],q[1],...q[n-1]
  m = 0;         //this->N_u*this->H_u;     // equality constraints g(x)
  nnz_jac_g = 0; // 2*m-1;           // nonzeros in the jacobian of the constraints g(x)
  index_style = C_STYLE;       // C index style for row/col entries

  return true;
}

bool NLMPC::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{ int i,j,k;

  for (j=0,k=0; j<N_u; j++, k++)
  { for (i=0; i<H_x; i++)   // lower bounds of the variables
    { x_l[k] = 0.0;
    }
  }

  for (j=0,k=0; j<N_u; j++, k++)
  { for (i=0; i<H_x; i++)   // upper bounds of the variables
    { x_u[k] = 13.0;
    }
  }
   
// bounds of the constraints
/*
  for (j=0,k=0; j<N_u; j++, k++)
  { for (i=0; i<H_x; i++)
    { g_l[k] = -2.0;
      g_u[k] = 2.0;
    }
  }
*/
  return true;
}


bool NLMPC::get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,Number* lambda)
{ int i, j, k;
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  for (j=0,k=0; j<N_u; j++, k++)
  { for (i=0; i<H_x; i++)
    { x[k] = this->u_applied(j);

    }
  }
  return true;
}

bool NLMPC::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) // cost function value at (x)
{ obj_value = this->costfct(n,(Number*)x);
  return true;
}

bool NLMPC::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) // cost function gradient value at (x)
{ Number xe[100];
  double eps=1e-6;
  int i, j, k;

  for (j=0,k=0; j<N_u; j++, k++)
  { for (i=0; i<H_u; i++)
    { xe[k]= x[k] +eps;
      grad_f[k] = (this->costfct(n,xe) -this->costfct(n,(Number*)x))/eps;
    }
  }
  return true;
}

bool NLMPC::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) // constraint value at (x)
{//return the value of the constraints: g(x).
  int i, j, k;

  // du max and min
/*
  for (j=0,k=0; j<N_u; j++, k++)
  { for (i=0; i<H_u; i++)   // lower bounds of the variables
    { if(k==0)
        g[k] = x[k]-u_applied(k);
      else if(k%H_u==0)
        g[k] = x[k]-u_applied(k);
      else
        g[k] = x[k]-x[i-1];
    }
  }
*/
  return true;
}

bool NLMPC::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values) //Return either the sparsity structure of the Jacobian of the constraints, or the values for the Jacobian of the constraints at the point $ x$ .
{ int i;
/*
  if (values == NULL)
  { // return the structure of the jacobian
    for (i=0; i<m; i++)   //setting the structure for the diagonal
    { iRow[i] = i;
      jCol[i] = i;
    }
    for (i=m; i<2*m-1; i++)   //setting the structure for the sub-diagonal
    { iRow[i] = i-m+1;
      jCol[i] = i-m;
    }
  }
  else
  { // return the values of the jacobian of the constraints
    for (i=0; i<m; i++)   //filling the values of the diagonal
      values[i] = 1.0;

    for (i=m; i<2*m-1; i++)
      values[i] = -1.0;
  }
*/
  return true;
}

bool NLMPC::eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values) //Return either the sparsity structure of the Hessian of the Lagrangian, or the values of the Hessian of the Lagrangian (9) for the given values for $ x$ , $ \sigma_f$ , and $ \lambda$ .
{ return true;
}

void NLMPC::finalize_solution(SolverReturn status,Index n, const Number* x,const Number* z_L,
                              const Number* z_U,Index m, const Number* g, const Number* lambda,
                              Number obj_value,
                              const IpoptData* ip_data,
                              IpoptCalculatedQuantities* ip_cq)
{ FILE * jFile;

  jFile = fopen( "ZReal_Input.txt", "a");//We write sent input to a txt file for plotting
  fprintf (jFile,"%f\n",x[0]);
  fclose(jFile);

  //Up to date of sent input for initialize next problem
  u_applied<< x[0], x[H_u];

  // MV prediction
  model_prediction(n,x); 

//read from real system
  calculate_output(n,x);
}