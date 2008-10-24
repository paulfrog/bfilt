// BFilt : A bayesian Filtering Library

//                     Copyright (C) 2008  Paul Frogerais

// The BFilt  Library is  free software: you  can redistribute  it and/or
// modify  it  under the  terms  of the  GNU  General  Public License  as
// published by  the Free  Software Foundation, either  version 3  of the
// License, or (at your option) any later version.

// This program  is distributed in the  hope that it will  be useful, but
// WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
// MERCHANTABILITY  or FITNESS  FOR  A PARTICULAR  PURPOSE.  See the  GNU
// General Public License for more details.

// You  should have received  a copy  of the  GNU General  Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

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
