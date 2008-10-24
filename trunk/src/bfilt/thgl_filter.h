///
/// @file   thgl_filter.h
/// @author paul <paul.frogerais@univ-rennes1.fr>
/// @date   Fri Sep 12 18:38:53 2008
/// 
/// @brief  The Taylor-Heun/Gauss-Legendre Filter (See Mazzoni 2007)
/// 
/// 
///

#ifndef THGL__FILTER__H
#define  THGL__FILTER__H

#include <bfilt/filter.h>


class THGL_Filter : public CD_Filter
{
public:

      THGL_Filter(void);
      THGL_Filter(Continuous_Discrete_Model *m);

protected :
      int  _update(const dcovector &Y);

};


#endif
