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
