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
