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
  \page page1 An AR process

  This is an example on the following auto regressive (AR) process : 
  \f[ X_k = 0.8 X_{k-1} + Wk    \f]
  where, \f$ X_{k} \in \mathcal{R} \f$, \f$ W_k \sim \mathcal{N}(0,0.1) \f$ \par

  This state is then observed by the output \f$ Y_k \in \mathcal{R} \f$ :
  \f[ Y_k = X_k + V_k \f]
  where \f$ V_k \sim \mathcal{N}(0,1) \f$
  
  \section sec1 Define the AR model

  First the ar process must be define as a sister class of Gaussian_Linear_Model.
  \include ar_process.h
  The Gaussian Linear Model are implemented in the following form :
  \f[ X_k = F X_{k-1} + f + G W_k \f]
  \f[ Y_k = H X_{k-1} + h + V_k \f]
  The constructor of AR_Process is then :
  \include ar_process.cpp

  \section sec2 The main program
  In the main program, the model will be first simulted with a specific 
  simulator for gaussian model (G_Simulator). Then the simulated output 
  sequence \f$ y_{0:N}\f$ is given to the input of a discrete-discrete 
  kalman filter (DD_Filter) to estimate the state \f$ \hat{X}_{0:k} \f$ \par
  First, all this objects are declared :
  \code
  int main(int argc, char **argv)
  {
      AR_Process model;            // The AR Model model

      G_Simulator sim(&model);      // The simulator
      
      DD_Kalman  filter(&model);   // The Kalman filter	

  \endcode

  Then 100 samples are simulated :

  \code 
  sim.Simulate(100);
  \endcode 

  The kalman filter is apply on the output sequence :
 
  \code 
  filter.Filtering(sim.Y);
  \endcode

  You can save the simulated sequences :
  \code 
      sim.Save_Y("output.dat");
      sim.Save_X("state.dat");
  \endcode
  and the estimated state :
  \code
      filter.Save_X("estimation.dat");
  \endcode
  With Gnuplot you can plot :
  \code
  plot 'state.dat' w l, 'estimation.dat' w l, 'output.dat'
  \endcode
  To obtain the following graph :
  \image html "ar_process.jpg" "Results" width=10
  \image latex ar_process.eps

  \section sec3 The CMakeList.txt
  \include CMakeLists.txt
*/

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
