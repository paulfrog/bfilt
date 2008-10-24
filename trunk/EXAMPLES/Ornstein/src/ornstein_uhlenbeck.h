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
