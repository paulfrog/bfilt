#include "van_der_pol.h"

Van_Der_Pol::Van_Der_Pol(void)
{
      lambda = 3.;

      Qw.resize(1);
      Qw(0,0)= 1.;
      Qv.resize(1);
      Qv(0,0)=0.1;
  
      X0.resize(2);
      X0(0) = 0.5;
      X0(1) = 0.5;
      R0.resize(2);
      R0.zero();
      R0(0,0)=0.;
      R0(1,1)=.1;
      Ts=.1;

}


dcovector Van_Der_Pol::Drift_Function(const dcovector & X)
{
      dcovector dX(X.l);

      dX(0) = X(1);
      dX(1) = lambda * (1. - X(0) * X(0)) * X(1) - X(0);

      return dX;
}

dgematrix Van_Der_Pol::J_Drift_Function(const dcovector & X)
{
      dgematrix F(X.l,X.l);

      F(0,0) = 0.;
      F(0,1) = 1.;
      F(1,0) = -2. * lambda * X(0) * X(1);
      F(1,1) = - lambda * X(0) * X(0);

      return F;
}
dcovector Van_Der_Pol::Observation_Function(const dcovector& X)
{
      dcovector Y(1);
      Y(0) = X(0);
      return Y;
}

dgematrix Van_Der_Pol::J_Observation_Function(const dcovector & X)
{
      dgematrix H(1,2);
      H(0,0) = 0.;
      H(0,1) = 1.;
      return H;
}
dgematrix Van_Der_Pol::Diffusion_Function(void)
{
      dgematrix G(2,1);
      G(0,0) = 0.;
      G(1,0) = 1.;

      return G;
}
