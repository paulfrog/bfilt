#ifndef __PLANE
#define __PLANE

#include <bfilt/gaussian_model.h>

class Plane : public Gaussian_Nonlinear_Model
{
      vector<double> Map;
      double xmin;
      double xmax;
      double ymin;
      double ymax;

      double sigv;
      double sigc;
public :
      Plane(const char *filename);

      dcovector State_Function(const dcovector &X, const dcovector &W);
      dcovector Observation_Function(const dcovector & X);
};


#endif
