#include <bfilt/simulator.h>
#include <bfilt/filter.h>
#include "ar_process.h"



int main(int argc, char **argv)
{
      AR_Process model;            // The AR model

      G_Simulator sim(&model);      // The simulator
      
      DD_Kalman  filter(&model);   // The Kalman filter	
      
      // Simulation 100 samples
      sim.Simulate(100);
        

      // Filtering from the simulated output sim.Y
      filter.Filtering(sim.Y);

      // Output Files for simulation
      sim.Save_Y("output.dat");
      sim.Save_X("state.dat");

      // Output File for filtering
      filter.Save_X("estimation.dat");
 
      return 0;
}
