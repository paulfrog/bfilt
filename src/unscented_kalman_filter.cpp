#include <bfilt/unscented_kalman_filter.h>

Unscented_Kalman_Filter::Unscented_Kalman_Filter(void)
{
}
Unscented_Kalman_Filter::Unscented_Kalman_Filter(Gaussian_Nonlinear_Model *m):GA_Filter::GA_Filter(m)
{
}


int Unscented_Kalman_Filter::_init(void)
{
      Gaussian_Nonlinear_Model *m=dynamic_cast<Gaussian_Nonlinear_Model *>(model);
      dcovector x0;
      dsymatrix r0;

      m->Get_Init_Parameters(x0,r0);
      M=x0;
      R=r0;

      if(sqrtm(m->Qw,sqrt_Qw))
            return 1;
      if(sqrtm(m->Qv,sqrt_Qv))
            return 1;
      int Nx=  M.l;
      int Nw= sqrt_Qw.n;
      int N= Nx+Nw;
      int i;
      lambda = 3. - N;

      w_0  = lambda / (N + lambda);
      w_0c = lambda / (N + lambda);// + (3. - alpha * alpha )
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


      U_Cov(sX,mx,sX,mx,Rp);
      Xp=mx;


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
      double phi = sqrt(N+lambda);

      if(sqrtm(R,sqrt_R))
            {
                  return 1;
            }
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
