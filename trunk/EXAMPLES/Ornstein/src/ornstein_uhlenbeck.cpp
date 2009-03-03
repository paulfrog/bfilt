#include "ornstein_uhlenbeck.h"

Ornstein_Uhlenbeck_Model::Ornstein_Uhlenbeck_Model(void)
{
      // parameters
      w = 4.;
      gamma = 2.;
      b = 8.;
      g = 2.;

      // Matrices of the state equation
      A.resize(2,2);
      A(0,0) = 0.;
      A(0,1) = 1.;
      A(1,0) = - (w*w);
      A(1,1) = -gamma;
      
      
      B.resize(2);
      B(0) = 0;
      B(1) = b;

      C.resize(2,1);
      C(0,0)=0.;
      C(1,0)=g;

      // Matrices of the Observation equation
      H.resize(1,2);
      H(0,0) = 1.;
      H(0,1) = 0.;

      h.resize(2);
      h.zero();

      Qw.resize(1);
      Qw.identity();
      Qw*=0.01;

      Qv.resize(1);
      Qv(0,0) = 0.001;

      // Sampling period
      Ts = 0.2;

      // Initial conditions
      R0.resize(2);
      R0.zero();
      R0(0,0)=0.;
      R0(1,1)=3.;
      
      X0.resize(2);
      X0.zero();
            

}
