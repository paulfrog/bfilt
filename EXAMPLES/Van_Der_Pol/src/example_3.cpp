#include <bfilt/simulator.h>
#include <bfilt/extended_kalman_filter.h>
#include "van_der_pol.h"



int main(int argc, char **argv)
{
      Van_Der_Pol model;              // The model

      CD_Simulator sim(&model);   // The simulator
      
      CD_Extended_Kalman_Filter  filter(&model,THGL);      // The filter
      
      // Simulation 40 seconds
      sim.Simulate(40.);
        

      // Filtering from the simulated output sim.Y
      filter.Filtering(sim.Y);

      // Output Files for simulation
      sim.Save_Y("output.dat");
      sim.Save_X("state.dat");

      // Output File for filtering
      filter.Save_X("estimation.dat");
 
      return 0;
}
