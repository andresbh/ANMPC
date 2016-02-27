//*** Main.c ***

#include <iostream>
#include "NMPCFormulation.hpp"
//#include "gnuplot_i.hpp"
//#include "EKF.h"

bool main(int argv, char* argc[])
{ /* 
  // plot ---------------------
  Gnuplot g1 = Gnuplot("lines");
  Gnuplot g2 = Gnuplot("lines");

  g1.set_style("steps");
  g2.set_style("steps");

  std::vector<double> meas1_v;
  std::vector<double> meas2_v;
  VectorXd yplot; 
  // plot ---------------------
  */
  VectorXi ublk;
  ublk.resize(2);
  ublk<<0,1;//,5,8,10;
  
  // Create an instance of a Non Linear MPC
  NMPC *mynmpc = new NMPC( 4, 4, 2, 2, 5, ublk, 1);
  // plot ------------
  //std::vector<double> ublk_v(ublk.data(), ublk.data()+ublk.size());
  //MatrixXd u_map = Map<MatrixXd>(mynmpc->u_applied.data(), mynmpc->N_u, mynmpc->H_u);
  //std::vector<double> u1_map(u_map.row(0).data(), u_map.row(0).data()+u_map.cols());
  //g1.plot_xy(ublk_v,u1_map,"u pred");
  //yplot.resize(mynmpc->N_y); // states;
  //yplot = Map<VectorXd,0,InnerStride<2>>((double*)mynmpc->measurements.data(),mynmpc->N_y,1);
  // plot -----------
  try
  { mynmpc->measurements(0)=8.2444;
    mynmpc->measurements(1)=19.0163;
    mynmpc->measurements(2)=4.3146;
    mynmpc->measurements(3)=8.8065;
    mynmpc->u_applied.fill(0.01);
    mynmpc->ref.row(0).fill(mynmpc->measurements(0));
    mynmpc->ref.row(1).fill(mynmpc->measurements(1));

    mynmpc->system_output(mynmpc->u_applied); // update from real plant
    int k=0;
	std::cout << "k ";
	std::cout << "r1 r2 ";
	std::cout << "m1 m2 ";
	std::cout << "x1 x2 x3 x4 ";
	std::cout << "xkf1 xkf2 xkf3 xkf4 xkf5 xkf6 ";
	//std::cout << "u1 u2 ";
	std::cout << "obj"<< std::endl;
	while(true)
    { k++;
      //Weight vectors
      mynmpc->Q.row(0).fill(1.0E5);  // system outputs
      mynmpc->Q.row(1).fill(1.0E5);

      mynmpc->R.fill(1.0E-2);        // system inputs
      mynmpc->L.row(0).fill(0.0); // system inputs linear term of cost function
      mynmpc->L.row(1).fill(0.0); // system inputs
      mynmpc->update();
      mynmpc->system_output(mynmpc->u_applied); // update from real plant

      if (!mynmpc->ok)
      { std::cout<<"\n\n*** Error during optimization!\n";
        system("PAUSE");
        return mynmpc->ok;
      }
      if(k> 500)
      { mynmpc->ref.row(0).fill(10);
        mynmpc->ref.row(1).fill(18);
      }
      std::cout<<k<<' ';
      std::cout<<mynmpc->ref.col(0).transpose() << ' ';
      std::cout<<mynmpc->measurements.transpose() << ' ';
      std::cout<<mynmpc->EKFS->x.transpose() << ' ';
      //std::cout<<mynmpc->solution->x.transpose() << ' ';
      std::cout<<mynmpc->solution->obj_value<<std::endl;

      // plot --------------
      //g1.reset_plot();
      //u_map = Map<MatrixXd>(mynmpc->u_applied.data(), mynmpc->N_u, mynmpc->H_u);
      //std::vector<double> u1_map(u_map.row(0).data(), u_map.row(0).data()+u_map.cols());
      //std::vector<double> u2_map(u_map.row(1).data(), u_map.row(1).data()+u_map.cols());
      //g1.plot_xy(ublk_v,u1_map,"u1 pred");
      //g1.plot_xy(ublk_v,u2_map,"u2 pred");
      /*
      yplot = Map<VectorXd,0,InnerStride<2>>((double*)mynmpc->measurements.data(),mynmpc->N_y,1);
      meas1_v.push_back(yplot(0));
      meas2_v.push_back(yplot(1));
      g1.reset_plot();
      g1.plot_x(meas1_v,"meas1");
      g2.reset_plot();
      g2.plot_x(meas2_v,"meas2");
      
      //g2.plot_x(meas2_v,"meas2");
      // plot --------------
	  */
    }
    
    delete mynmpc;

    system("PAUSE");
    // As the SmartPtrs go out of scope, the reference count
    // will be decremented and the objects will automatically 
    // be deleted.
  }
  catch(std::exception const& e)
  { std::cout << "There was an error: " << e.what() << std::endl;
    system("PAUSE");
  }
  return mynmpc->ok;
  
}