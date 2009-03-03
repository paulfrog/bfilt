#include <bfilt/local_linearization_filter.h>


LL_Filter::LL_Filter(void)
{
      Likelihood=0.;
}

LL_Filter::LL_Filter(Continuous_Discrete_Model *m):CD_Filter::CD_Filter(m)
{
      this->model=m;
      Init();

}

int  LL_Filter::_update(const dcovector &Y)
{
      Continuous_Discrete_Model *m=dynamic_cast<Continuous_Discrete_Model *>(model);
      const double _2_PI_= 6.283185307;

      double delta = m->Ts;

      dgematrix J=m->J_Drift_Function(M);
      dcovector f=m->Drift_Function(M);
      dgematrix J1=i(J);
      dgematrix F = expm(J*delta);
      dgematrix G = m->Diffusion_Function();
      dgematrix H;
      dgematrix Qi;
      dcovector I;
      dgematrix K;
      dsymatrix Qw=m->Qw*delta;
      dsymatrix Qv=m->Qv;


      Xp = M + (J1*F - J1)*f;

      H=m->J_Observation_Function(Xp);  
      Rp = F * R * t(F) + G * (Qw * t(G));
  
      Qi =  H * (Rp * t(H)) + Qv;

      K = Rp* t(H) * i(Qi);
      I= (Y - m->Observation_Function(Xp));
      M= Xp + K * I;
      R= Rp - K * H * Rp;

      Likelihood += - 0.5 * I.l * log(_2_PI_) - 0.5 * log(det(Qi)) - 0.5 * CPPL::t(I) * CPPL::i(Qi) * I;
      return 0;
}

