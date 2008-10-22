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

#include <bfilt/extended_kalman_filter.h>

Extended_Kalman_Filter::Extended_Kalman_Filter(void)
{
      Likelihood=0.;
}


Extended_Kalman_Filter::Extended_Kalman_Filter(Gaussian_Nonlinear_Model  *m):
      GA_Filter::GA_Filter(m)
{
}


int Extended_Kalman_Filter::_update(const dcovector &Y)
{
      Gaussian_Nonlinear_Model *m=dynamic_cast<Gaussian_Nonlinear_Model *>(model);
 
      const double _2_PI_= 6.283185307;

      dsymatrix Qw=m->Qw;
      dsymatrix Qv=m->Qv;
      dcovector W(Qw.n);
      W.zero();
      dgematrix F;
      dgematrix G;
      dgematrix H;
      dgematrix Qi;
      dcovector I;
      dgematrix K;

      m->Get_Linear_Parameters(M,W,F,G,Xp);
      

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




CD_Extended_Kalman_Filter::CD_Extended_Kalman_Filter(void)
{
}

CD_Extended_Kalman_Filter::CD_Extended_Kalman_Filter(Continuous_Discrete_Model *m, const int & sh)
{
      model=m;
      Scheme=sh;
}
int  CD_Extended_Kalman_Filter::_update(const dcovector &Y)
{
      Continuous_Discrete_Model *m=dynamic_cast<Continuous_Discrete_Model *>(model);

      const double _2_PI_= 6.283185307;
      
      dgematrix H;  // observation matrix
      dgematrix Qi;
      dcovector I;  // innovation
      dgematrix K;  // kalman gain

      dsymatrix Qv=m->Qv;

      switch(Scheme)
            {
            case EULER :
                  _euler_prediction(Xp,Rp);
                  break;
            case THGL :
                  _thgl__prediction(Xp,Rp);
                  break;
            case RK4 :
                  _rk4___prediction(Xp,Rp);
                  break;
            case HEUN :
                  _heun__prediction(Xp,Rp);
            default :
                  break;
            }
      
            
      
  
      H=m->J_Observation_Function(Xp);  

      Qi =  H * (Rp * t(H)) + Qv;

      K = Rp* t(H) * i(Qi);
      I= (Y - m->Observation_Function(Xp));
      M= Xp + K * I;
      R= Rp - K * H * Rp;

      Likelihood += - 0.5 * I.l * log(_2_PI_) - 0.5 * log(det(Qi)) - 0.5 * CPPL::t(I) * CPPL::i(Qi) * I;
      
      return 0;
}

void CD_Extended_Kalman_Filter::_thgl__prediction(dcovector &Xp, dgematrix &p)  // Prediction ref : T. Mazzoni (2007)
{
      Continuous_Discrete_Model *m=dynamic_cast<Continuous_Discrete_Model *>(model);
      dgematrix Id(M.l); // identity matrix
      dcovector f;  // drift function in X
      dgematrix J;  // Jacobian in X
      dgematrix J1; // Jacobian in X1
      dcovector X1; // intermeditate prediction at t+delta/2
      dgematrix G ; // Diffusion term
      double delta = m->Ts;
      dsymatrix Qw=m->Qw*delta;
      dgematrix A;
      dgematrix B;


      Id.identity();
      f=m->Drift_Function(M);
      J=m->J_Drift_Function(M);
      G= m->Diffusion_Function();
      Xp = M + CPPL::i(Id - J*(delta/2.)) *f * delta;
      X1 = 0.5 * (M + Xp - J * f * delta*delta / 4.);
      J1=m->J_Drift_Function(X1);
      A = CPPL::i(Id - J1 * delta / 2.);
      B= A*(Id + J1 * delta / 2.);
      p = A * R * t(A) + B *(G * (Qw * t(G)))*t(B);
}

void CD_Extended_Kalman_Filter::_euler_prediction(dcovector &Xp, dgematrix &p)
{
      Continuous_Discrete_Model *m=dynamic_cast<Continuous_Discrete_Model *>(model);
      double delta = m->Ts;
      dgematrix J1;
      dcovector dx1;
      dgematrix K1;

      dgematrix G = m->Diffusion_Function();

      dgematrix Q = G*(m->Qw*t(G));


      dx1 =  m->Drift_Function(M);
      Xp= M + delta * ( dx1 );

      J1 = m->J_Drift_Function(M);
      K1 = J1 * R + R * t(J1) + Q;
      p = R + delta * K1;     
}

void CD_Extended_Kalman_Filter::_rk4___prediction(dcovector &Xp, dgematrix &p)
{
      Continuous_Discrete_Model *m=dynamic_cast<Continuous_Discrete_Model *>(model);
      dgematrix G = m->Diffusion_Function();
      dgematrix Q=G*(m->Qw*t(G));
      double delta = m->Ts;

      dcovector  dx1, dx2, dx3, dx4;

      dgematrix K1, K2, K3, K4;

      dgematrix J1,J2,J3,J4;

      dx1 =  m->Drift_Function(M);

      dx2 =  m->Drift_Function(M + 0.5*delta *dx1);

      dx3 =  m->Drift_Function(M + 0.5*delta *dx2);

      dx4 =  m->Drift_Function(M + delta *dx3); 

      Xp= M + delta * (dx1+2.*dx2+2.*dx3+dx4)/6.;

      J1 = m->J_Drift_Function(M);
      J2 = m->J_Drift_Function(M + 0.5 * delta * dx1);
      J3 = m->J_Drift_Function(M + 0.5 * delta * dx2);
      J4 = m->J_Drift_Function(M + delta * dx3);

      K1 = J1 * R + R * t(J1) + Q;
      K2 = J2 * (R + 0.5*delta * K1) + (R + 0.5*delta * K1) * t(J2) + Q;
      K3 = J3 * (R + 0.5*delta * K2) + (R + 0.5*delta * K2) * t(J3) + Q;
      K4 = J4 * (R + delta * K3) + (R + delta * K3) * t(J4) + Q;

      p = R + delta *(K1 + 2.*K2 + 2.*K3 + K4) / 6.;
}

void CD_Extended_Kalman_Filter::_heun__prediction(dcovector &Xp, dgematrix &p)
{

      Continuous_Discrete_Model *m=dynamic_cast<Continuous_Discrete_Model *>(model);
      dgematrix G = m->Diffusion_Function();
      dgematrix Q=G*(m->Qw*t(G));
      double delta = m->Ts;

      dcovector dx1, dx2;

      dgematrix K1, K2;

      dgematrix J1, J2;

      dx1 =  m->Drift_Function(M);


      dx2 =  m->Drift_Function(M + delta *dx1); 

      Xp= M + delta * (dx1+dx2)/2.;

      J1 = m->J_Drift_Function(M);
      J2 = m->J_Drift_Function(M + delta * dx1);

      K1 = J1 * R + R * t(J1) + Q;
      K2 = J2 * (R + delta * K1) + (R + delta * K1) * t(J2) + Q;

      p = R + delta * (K1 + K2) /2.;
}
