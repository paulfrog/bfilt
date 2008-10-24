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

#ifndef __ORNSTEIN_UHL__
#define __ORNSTEIN_UHL__

#include <bfilt/gaussian_model.h>

class Ornstein_Uhlenbeck_Model : public Linear_CD_Model
{
public :
      double gamma;
      double w;
      double b;
      double g;
      Ornstein_Uhlenbeck_Model(void);


};

#endif
