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

#include <bfilt/unscented_kalman_filter.h>

Unscented_Kalman_Filter::Unscented_Kalman_Filter(Gaussian_Nonlinear_Model *m):GA_Filter::GA_Filter(m)
{

      model=m;
      alpha = 0.5;
}

Unscented_Kalman_Filter::Unscented_Kalman_Filter(Gaussian_Nonlinear_Model *m, const double & a):GA_Filter::GA_Filter(m)
{

      model=m;
      alpha = a;
}

int Unscented_Kalman_Filter::_init(void)
{
      Gaussian_Nonlinear_Model *m=dynamic_cast<Gaussian_Nonlinear_Model *>(model);
      dcovector x0;
      dsymatrix r0;
      m->Get_Init_Parameters(x0,r0);
      M=x0;
      R=r0;
      if(cholesky(m->Qw,sqrt_Qw))
            return 1;
      if(cholesky(m->Qv,sqrt_Qv))
            return 1;
      int Nx=  M.l;
      int Nw= sqrt_Qw.n;
      int N= Nx+Nw;
      int i;
      double lambda = alpha * alpha * N - N;

      w_0  = lambda / (N + lambda);
      w_0c = lambda / (N + lambda) + (3. - alpha * alpha );
      w    = 1. /(2.*(N+ lambda));

      sX.resize(2*N+1);
      sW.resize(2*N+1);
      sY.resize(2*N+1);

      Likelihood = 0.;
      return 0;
}


int Unscented_Kalman_Filter::_update(const dcovector &Y)
{
      Gaussian_Nonlinear_Model *m=dynamic_cast<Gaussian_Nonlinear_Model *>(model);

      const double _2_PI_= 6.283185307;

      int N,i;
      dgematrix cyy,cxy,k;
      dcovector my,mx;
      dcovector I;
      if(SP_Init())
            return 1;

      N=sX.size();

      for(i=0;i<N;i++)
	    {
		  sX[i]=m->State_Function(sX[i],sW[i]);
	    }

      U_Mean(sX,mx);
      Xp=mx;

      U_Cov(sX,mx,sX,mx,Rp);

      for(i=0;i<N;i++)
	    {
		  sY[i]=m->Observation_Function(sX[i]);
	    }

      U_Mean(sY,my);
      U_Cov(sY,my,sY,my,cyy);

      cyy=cyy + m->Qv;

      U_Cov(sX,mx,sY,my,cxy);

      k=cxy*CPPL::i(cyy);

      I = (Y-my);
      M = mx + k * I;

      R = Rp - k * cyy * t(k);


      Likelihood += -0.5* I.l * log(_2_PI_) - 0.5*log(det(cyy)) - 0.5 * CPPL::t(I) * CPPL::i(cyy) * I;

      return 0;
}

int Unscented_Kalman_Filter::SP_Init(void)
{
      int Nx=  M.l;
      int Nw= sqrt_Qw.n;
      int N= Nx+Nw;
      int i;
      dgematrix  sqrt_R;
      double phi = alpha*sqrt(N);

      if(cholesky(R,sqrt_R))
            return 1;

      sX[0]=M;
      sW[0].zero();

      for(i=1;i<1+Nx;i++){
	    sX[i]=M + phi*get_column(sqrt_R,i-1);
	    sW[i].zero();
      }
      for(i=Nx+1;i<1+Nx+Nw;i++){
	    sX[i]=M;
	    sW[i]=phi*get_column(sqrt_Qw,i-(Nx+1));
      }
      for(i=1+Nx+Nw;i<1+2*Nx+Nw;i++){
	    sX[i]=M-phi*get_column(sqrt_R,i-(1+Nx+Nw));
	    sW[i].zero();
      }
      for(i=1+2*Nx+Nw;i<1+2*Nx+2*Nw;i++){
	    sX[i]=M;
	    sW[i]=-phi*get_column(sqrt_Qw,i-(1+2*Nx+Nw));
      }
      return 0;
}

int Unscented_Kalman_Filter::U_Mean(const vector<dcovector> &sP, 
			    dcovector & mean)
{
      int N=sP[0].l;
      int i;
      mean.resize(N);
      mean=sP[0]*w_0;
      N=sP.size();
      for (i=1;i<N;i++)
	    mean +=w * (sP[i]);
}

int Unscented_Kalman_Filter::U_Cov(const vector<dcovector> &sP1, 
			   const dcovector &m1,
			   const vector<dcovector> &sP2, 
			   const dcovector &m2,
			   dgematrix & cov)
{
      int N=sP1[0].l;
      int M=sP2[0].l;
      int i;
      cov.resize(N,M);
      cov.zero();
      N=sP1.size();
      drovector Xt = CPPL::t(sP2[0] - m2);
      dcovector X = (sP1[0] - m1);
      cov =  X* Xt;
      cov*=w_0c;
      for(i=1;i<N;i++)
	    cov += w * ((sP1[i] - m1) * CPPL::t(sP2[i] - m2));
}

dcovector get_column(const dgematrix &M, const int & k)
{
      dcovector X;
      int i;
      X.resize(M.m);
      for (i=0;i<M.m;i++)
	    X(i)=M(i,k);
      return X;
}
