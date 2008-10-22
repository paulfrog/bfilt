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

///
/// @file   local_linearization_filter.h
/// @author paul <paul.frogerais@univ-rennes1.fr>
/// @date   Fri Sep 12 18:37:04 2008
/// 
/// @brief  Implementation of the local linearization filter (See Ozaki 1993)
/// 
/// 
///


#ifndef LL__FILTER__H
#define  LL__FILTER__H

#include <bfilt/filter.h>

class LL_Filter : public CD_Filter
{
public:

      LL_Filter(void);
      LL_Filter(Continuous_Discrete_Model *m);

protected :
      int  _update(const dcovector &Y);
};


#endif
