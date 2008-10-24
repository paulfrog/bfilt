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

#ifndef __VAN_DER_POL
#define __VAN_DER_POL


#include <bfilt/gaussian_model.h>


class Van_Der_Pol : public Continuous_Discrete_Model
{
public :
      double lambda;
      
      Van_Der_Pol(void);
      dcovector Drift_Function(const dcovector & X);
      dgematrix J_Drift_Function(const dcovector & X);
      dcovector Observation_Function(const dcovector& X);
      dgematrix J_Observation_Function(const dcovector & X);
      dgematrix Diffusion_Function(void);
};


#endif
