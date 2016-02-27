#include "Cost.hpp"

COSTF::COSTF()
{
}

COSTF::~COSTF()
{
}

COSTF::COSTF(const int N_x0, const int N_y0, const int N_z0, const int N_u0, const int H_x0, const int H_u0, const VectorXi ublk0 , const double h0)
 : model(N_y0, N_z0, N_u0, true, false)
{ //int n_thread = NUMBER_THREADS;   // number of threads in parallel regions
	//omp_set_dynamic(0);              // off dynamic thread adjust
	//omp_set_num_threads(n_thread);   // set the number of threads 

  N_x = N_x0;
  N_y = N_y0;
  N_z = N_z0;
  N_u = N_u0;
  H_x = H_x0;
  H_u = H_u0;
  MatrixXd C0;
  VectorXd P0;
  VectorXd Q0;
  VectorXd R0;

  ublk = ublk0;

  // Initial states
  measurements.resize(N_y+N_z,1);
  measurements(0)=8;
  measurements(1)=19;
  measurements(2)=4;
  measurements(3)=9;

  ADvector y3; y3.resize(N_y+N_z);
  ADvector u3; u3.resize(N_u);

  C0.resize(N_y+N_z,N_y+N_z);
  C0 << MatrixXd::Identity(N_x+N_z, N_x+N_z);
  // y1 controlled // y2 controlled
  P0.resize(N_x+N_z);
  //  x1  x2  x3  x4  r1  r2
  P0<<1,1,1,1,1,1;
  P0*=1e-2;

  Q0.resize(N_x+N_z);
  Q0 << 1.0E-2,1.0E-2,1.0E0,1.0E0,1.0E0,1.0E0;
  Q0*=1e0;

  R0.resize(N_x+N_z);
  R0 << 1.0E-5,1.0E-5,1.0E0,1.0E0,1.0E-1,1.0E-1;
  R0*=1e-2;
  //R0 << 1.0E0,1.0E0,1.0E0,1.0E0,1.0E0,1.0E0;

  EKFS = new EKF<COSTF, AD<double>, AD<double>>
         ( measurements.cast<AD<double>>(),
           y3, u3, h0, 
           C0.cast<AD<double>>(),
           N_x, N_z,
           P0.cast<AD<double>>().asDiagonal(), 
           Q0.cast<AD<double>>().asDiagonal(), 
           R0.cast<AD<double>>().asDiagonal()
         );

  this->ref.resize(N_z,H_x);
  this->Q.resize(N_z,H_x);
  this->R.resize(N_u,H_u);
  this->L.resize(N_u,H_u);
  this->x_KF.resize(N_x+N_z,H_x+1);
}

AD<double> COSTF::costfct(const ADvector u)
{ int i, k; // counters
  AD<double> costvalue = 0.0;

  Matrix<ADdouble, Dynamic, Dynamic> u0_ = Map<Matrix<ADdouble, Dynamic, Dynamic>>((ADdouble*)u.data(), N_u, H_u);
  Matrix<ADdouble, Dynamic, 1> uapplied0_ = this->u_applied.cast<ADdouble>().head(N_u);
  this->x_KF.col(0) = EKFS->x;
  //EKFS->update(this->measurements, u0_.col(0));
  
  //std::cout<<std::endl<<"x_kf1: "<<EKFS->x.transpose()<<std::endl;
  //std::cout<<"meas: "<<this->measurements.transpose()<<std::endl;
  //std::cout<<"kf: "<<EKFS->y.transpose()<<std::endl;
  //std::cout<<"u0: "<<u0_.col(0).transpose()<<std::endl;
  costvalue +=(this->measurements.head(N_z)-this->ref.col(0)).transpose()*Q.col(0).asDiagonal()*(this->measurements.head(N_z)-this->ref.col(0));
  //# pragma omp parallel for
  for (i=0,k=0; i<this->H_x; i++)
  { this->x_KF.col(i+1) = EKFS->update(this->x_KF.col(i), EKFS->y, u0_.col(k));
    //std::cout<<i<<": ref:  "<<this->ref.col(i)<<std::endl;
    // Call model to get derivatives evaluated at these values.
    //this->x_KF.segment(N_x, N_z) = ref.col(i);
    //std::cout<<i<<": x_KF:  "<<this->x_KF.col(i+1).transpose()<<std::endl;
    //std::cout<<i<<": u0:  "<<u0_.col(k)<<std::endl;
    //std::cout<<i<<": x_KF:  "<<this->x_KF.transpose()<<std::endl;

    //std::cout<<i<<": u0:  "<<u0_.col(k).transpose()<<std::endl;
    //std::cout<<i<<": uapplied0_:  "<<uapplied0_.transpose()<<std::endl;
    costvalue += (u0_.col(k)-uapplied0_).transpose()*(R.col(k).asDiagonal()*(u0_.col(k)-uapplied0_)); 
    //std::cout<<"cost1:" <<u0_.col(k).transpose()*(R.col(k).asDiagonal()*u0_.col(k))<<std::endl;
    costvalue += (L.col(k).transpose()*u0_.col(k)); 
    //std::cout<<"cost2:" <<(L.col(k).transpose()*u0_.col(k))<<std::endl;
    //costvalue +=(EKFS->x.head(N_z)).transpose()*Q.col(i).asDiagonal()*(EKFS->x.head(N_z));
    //std::cout<<"y_c:" <<EKFS->y.head(N_z).transpose()<<" ref: "<< ref.col(i).transpose()<<" u: "<<u0_.col(k).transpose()<<std::endl;
    //std::cout<<"cost3:" <<((EKFS->y.head(N_z)-this->ref.col(i)).transpose()*Q.col(i).asDiagonal()*(EKFS->y.head(N_z)-this->ref.col(i))).transpose()<<std::endl;
    //costvalue +=(EKFS->y.head(N_z)-this->ref.col(i)).transpose()*Q.col(i).asDiagonal()*(EKFS->y.head(N_z)-this->ref.col(i));
    costvalue +=(this->x_KF.col(i+1).head(N_z)-this->ref.col(i)).transpose()*Q.col(i).asDiagonal()*(this->x_KF.col(i+1).head(N_z)-this->ref.col(i));
    
    if((ublk(k)==i)&&(k<ublk.size()-1))
      k++;
  }
  EKFS->x = this->x_KF.col(0);
  //std::cout<<std::endl<<"x_kf: "<<std::endl<<this->x_KF<<std::endl;
  
  //EKFS->update(this->measurements, u0_.col(0));

  //std::cout<<"obj: "<<costvalue<<std::endl;
  return costvalue;
}

int COSTF::system_output(Dvector u)
{ Dvector y; y.resize(this->N_y); // states;
  Dvector dy;dy.resize(this->N_y);// derivatives;
  y << Map<Dvector,0,InnerStride<2>>((double*)measurements.data(),N_y,1);
  MatrixXd u_map = Map<MatrixXd>(u.data(), this->N_u, this->H_u);

  // Update state with Euler time discr
  // Simulate using simple forward Differences
  //std::cout<<"u_mod: "<<u_map.col(0)<<std::endl;
  //std::cout<<"y_out: "<<y.transpose()<<std::endl;
  realModel( dy.data(), u_map.col(0).data(), y.data()); // apply only first control
  // Call model to get derivatives evaluated at these values.
  y += EKFS->h*dy;
  // update system outputs
  this->measurements.head(N_y)<<y.cast<ADdouble>();
  this->measurements.segment(N_y,N_z)<<this->ref.col(0); 

  return 0;
}

int COSTF::model_prediction(const Dvector u)
{ int i, k;    // Counters
  VectorXd x = VectorXd::Zero(N_x) ;  // states
  VectorXd dx = VectorXd::Zero(N_x); // derivatives
  
  // Update state with Euler time discr
  for (i=0, k=0; i<H_x; i++)
  { // Simulate using simple forward Differences
    //EKFS->update(this->measurements, u0_.col(k));
    
    // Call model to get derivatives evaluated at these values.
    x += 0.1*dx;
    if(ublk(k)==i)
      k++;
  }

  return 0;
}


int COSTF::operator()(const Matrix<AD<double>, Dynamic, 1> x, 
                      const Matrix<AD<double>, Dynamic, 1> u, 
                      const Matrix<AD<double>, Dynamic, 1> d, 
                            Matrix<AD<double>, Dynamic, 1> &fvec)
{ modelcalc(x, u, d, fvec );
  return 1;
}