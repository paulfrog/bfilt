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

/*! 
  \page page1 A linear Discrete-Discrete Exemple

  This is a little exemple on the following auto regressive (AR) process : 
  \f[ X_k = 0.8 X_{k-1} + Wk    \f]
  where, \f$ X_{k} \in \mathcal{R} \f$, \f$ W_k \sim \mathcal{N}(0,1) \f$ \par

  This state is then observed by the output \f$ Y_k \in \mathcal{R} \f$ :
  \f[ Y_k = X_k + V_k \f]
  where \f$ V_k \sim \mathcal{N}(0,0.001) \f$
  
  \section sec1 Define the AR model

  First the ar process must be define as a sister class of Gaussian_Linear_Model.
  \include ar_process.h
  The Gaussian Linear Model are implemented in the following form :
  \f[ X_k = F X_{k-1} + f + G W_k \f]
  \f[ Y_k = H X_{k-1} + h + V_k \f]
  The constructor of AR_Process is then :
  \include ar_process.cpp

  \section sec2 The main program
  \subsection subsec1 Simulation of the process
  \subsection subsec2 Kalman Filtering of the process
*/
#include <bfilt/simulator.h>
#include <bfilt/filter.h>
#include "ar_process.h"



int main(int argc, char **argv)
{
      AR_Process model;            // A nonlinear continuous discrete (CD) model

      Simulator sim(&model);      // A simulator
      
      THGL_Filter  filter(&model);   // A CD EKF filter	
      
      // Simulation 100 samples
      sim.Simulate(100);
        
      // Output Files for simulation
      sim.Save_Y("output.dat");
      sim.Save_X("state.dat");

      // Filtering from the simulated output sim.Y
      filter.Filtering(sim.Y);

      // Output File for filtering
      filter.Save_X("estimation.dat");
 
      return 0;
}
