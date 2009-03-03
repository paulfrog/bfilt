#include <bfilt/simulator.h>
#include <bfilt/filter.h>
#include "ornstein_uhlenbeck.h"



int main(int argc, char **argv)
{
      Ornstein_Uhlenbeck_Model model; // The model

      LTI_CD_Simulator sim(&model);   // The simulator
      
      CD_Kalman  filter(&model);      // The Kalman filter	
      
      // Simulation 10 seconds
      sim.Simulate(10.);
        

      // Filtering from the simulated output sim.Y
      filter.Filtering(sim.Y);

      // Output Files for simulation
      sim.Save_Y("output.dat");
      sim.Save_X("state.dat");

      // Output File for filtering
      filter.Save_X("estimation.dat");
 
      return 0;
}
