#include "ar_process.h"

AR_Process::AR_Process(void)
{
      // State Equation
      F.resize(1,1);
      F(0,0) = 0.8;

      f.resize(1);
      f(0) = 0.;
      
      G.resize(1,1);
      G.identity();
      
      Qw.resize(1);
      Qw(0,0)= 0.1;

      // Observation noise
      H.resize(1,1);
      H(0,0) = 1;

      h.resize(1);
      h(0) = 0.;

      Qv.resize(1);
      Qv(0,0)=1;
  
      // Init state 
      X0.resize(1);
      X0(0) = 10.;

      R0.resize(1);
      R0.zero();
}
