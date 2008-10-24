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
