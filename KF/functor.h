///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic, int N = 0>
struct Functor
{ typedef _Scalar Scalar;
  enum {
      InputsAtCompileTime = NX,
      ValuesAtCompileTime = NY
  };
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  //typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType; 

  int m_inputs, m_values, m_Nl;
  double m_l;
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> m_M, m_B, m_Y, m_E, m_D;

  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime), m_Nl(Nl) {}
  Functor(int inputs, int values, double ld, int Nl) : m_inputs(inputs), m_values(values), m_l(ld), m_Nl(Nl) {
    const Scalar alpha = 1; // numerical stability choose between (0-1)
    m_M.Zero(m_Nl, m_Nl);
    m_M(0, 0) =(m_l/3 -m_l*alpha/4);
    m_M(m_Nl-1, m_Nl-1) =m_l/3 +m_l*alpha/4;
    m_M(0, 1) =m_l/6 -m_l*alpha/4;
    m_M(m_Nl-1, (m_Nl-2)) =m_l/6 +m_l*alpha/4;
    if (m_Nl > 2)
    {  for (int i =1; i <m_Nl-2; i++)
       { m_M(i, i -1) =m_l/6 +m_l*alpha/4;
         m_M(i, i   ) =2*m_l/3;
         m_M(i, i +1) =m_l/6 -m_l*alpha/4;
       }
    }
    m_B.Zero(m_Nl, m_Nl);
  }

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

};

struct my_functor : Functor<double>
{ my_functor(void): Functor<double>(2,2,10.,5){}
  
  void operator()( const Matrix<double, Dynamic, 1> x, 
                   const Matrix<double, Dynamic, 1> u, 
                         Matrix<double, Dynamic, 1> fvec)
  { // Implement y = 10*(x0+3)^2 + (x1-5)^2
    fvec(0) = 10.0*pow(x(0)+3.0, 2.0) +  pow(x(1)-5.0, 2.0);
    fvec(1) = 2.0*u(0);

    return;
  }
};

///////////////////////////////////////////////////////////////////////////////