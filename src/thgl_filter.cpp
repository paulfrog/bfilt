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

#include <bfilt/thgl_filter.h>


THGL_Filter::THGL_Filter(void)
{
}

THGL_Filter::THGL_Filter(Continuous_Discrete_Model *m):CD_Filter::CD_Filter(m)
{
      this->model=m;
      Init();

}

int  THGL_Filter::_update(const dcovector &Y)
{
      Continuous_Discrete_Model *m=dynamic_cast<Continuous_Discrete_Model *>(model);

      const double _2_PI_= 6.283185307;

      double delta = m->Ts;
      dgematrix A;
      dgematrix B;
      dgematrix J;  // Jacobian in X
      dgematrix J1; // Jacobian in X1
      dcovector f;  // drift function in X
      dcovector X1; // intermeditate prediction at t+delta/2
      dgematrix G ; // Diffusion term
      dgematrix H;  // observation matrix
      dgematrix Qi;
      dcovector I;  // innovation
      dgematrix K;  // kalman gain
      dsymatrix Qw=m->Qw*delta;
      dsymatrix Qv=m->Qv;
      dgematrix Id(M.l); // identity matrix

      Id.identity();
      f=m->Drift_Function(M);
      // Prediction ref : T. Mazzoni (2007)
      J=m->J_Drift_Function(M);
      G= m->Diffusion_Function();
      Xp = M + CPPL::i(Id - J*(delta/2.)) *f * delta;
      X1 = 0.5 * (M + Xp - J * f * delta*delta / 4.);
      J1=m->J_Drift_Function(X1);
      A = CPPL::i(Id - J1 * delta / 2.);
      B= A*(Id + J1 * delta / 2.);
      
            
      Rp = A * R * t(A) + B *(G * (Qw * t(G)))*t(B);
  
      H=m->J_Observation_Function(Xp);  

      Qi =  H * (Rp * t(H)) + Qv;

      K = Rp* t(H) * i(Qi);
      I= (Y - m->Observation_Function(Xp));
      M= Xp + K * I;
      R= Rp - K * H * Rp;

      Likelihood += - 0.5 * I.l * log(_2_PI_) - 0.5 * log(det(Qi)) - 0.5 * CPPL::t(I) * CPPL::i(Qi) * I;
      return 0;
}
