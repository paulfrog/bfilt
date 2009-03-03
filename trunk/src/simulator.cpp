#include <bfilt/simulator.h>

Simulator::Simulator(void)
{
      gsl_rng_env_setup(); 
      r = gsl_rng_alloc (gsl_rng_taus);
      gsl_rng_set(r,time(0)); 
      b=NULL;
      model=NULL;
}

Simulator::~Simulator(void){
      gsl_rng_free(r);
}

void Simulator::Set_Seed(const int & s)
{
      gsl_rng_set(r,s); 
}
void  Simulator::_update(void)
{
      int N=X.size();
      if (N ==0)
            {
                  if(model)
                        model->Clear();
                  X.push_back(Draw_Init());
                  Y.push_back(Draw_Observation(X[0]));
            }
      else
            {
                  if(model)
                        model->Update();
                  X.push_back(Draw_Transition(X[N-1]));
                  Y.push_back(Draw_Observation(X[N-1]));
            }
}

void  Simulator::Update(void)
{
      _update();
}

void Simulator::Simulate(const int & N)
{
      
      X.resize(N);
      Y.resize(N);
      X[0]=Draw_Init();
      Y[0]=Draw_Observation(X[0]);
      if(model)
            model->Clear();

      for (int i=1;i<N;i++){
            if(model)
                  model->Update();
	    X[i]=Draw_Transition(X[i-1]);
	    Y[i]=Draw_Observation(X[i]);
      }
}
void Simulator::Clear(void)
{
      Y.clear();
      X.clear();
}

int Simulator::Save_Y(const char *filename)
{
      ofstream file(filename);
      int i,j;
      int N=Y.size();

      if(file){
	    for (i=0;i<N;i++)
                  {
                        for (j=0; j<Y[0].l; j++)
                              file<<Y[i](j)<<" ";
                        file<<endl;
	    }
	    file.close();
      }
      
      else
	    {
		  cout<<"ERROR : Pb to access file"<<endl;
		  return 1;
	    }
      return 0;

}
int Simulator::Save_X(const char *filename)
{
      ofstream file(filename);
      int i,j;
      int N=X.size();
      
      if(file){
	    for (i=0;i<N;i++)
                  {
                        for (j=0; j<X[0].l; j++)
                              file<<X[i](j)<<" ";
                        file<<endl;
	    }
	    file.close();
      }
      
      else
	    {
		  cout<<"ERROR : Pb to access file"<<endl;
		  return 1;
	    }
      return 0;

}


CD_Simulator::CD_Simulator(void)
{
}

void CD_Simulator::Set_Alpha(const int & al)
{
      _a = al;
      Dx =  Dy / _a; 
}

CD_Simulator::CD_Simulator(Continuous_Discrete_Model *cd_m, const int & s, const int &alpha)
{
      _a= alpha;
      Dy=cd_m->Ts;
      Dx =  Dy / _a; 
      scheme = s;
      model = cd_m;
}


dcovector CD_Simulator::draw_state(const dcovector &X)
{
      Continuous_Discrete_Model *cd_m=dynamic_cast<Continuous_Discrete_Model *>(model);
      
      dgematrix Q= cd_m->Qw * Dx;
      dcovector w;
      dcovector U(X.l), dx1(X.l), dx2(X.l), dx3(X.l), dx4(X.l);
      dgematrix G=cd_m->Diffusion_Function();
      dgematrix J;
      dcovector f;
      dgematrix J1;
      dgematrix R;
      if(b)
            w=b(this,r);
      else
            multivariate_normal_draw(w,r,Q);

      switch(scheme)
            {

            case EULER :
                  U= X + Dx*cd_m->Drift_Function(X) +G*w;
                  break;

            case HEUN  :
                  dx1 = cd_m->Drift_Function(X) + G*w/Dx;

                  dx2 = cd_m->Drift_Function(X + dx1*Dx) + G*w/Dx; 
  
                  U=X+ (dx1+dx2)*Dx/2.;
                  break;

            case SRK4  :
                  dx1 =  cd_m->Drift_Function(X) + G*w/Dx;

                  dx2 =  cd_m->Drift_Function(X + 0.5*Dx *dx1)+G*w/Dx;

                  dx3 =  cd_m->Drift_Function(X + 0.5*Dx *dx2)+G*w/Dx;

                  dx4 =  cd_m->Drift_Function(X + Dx *dx3) + G*w/Dx; 
  
                  U=X+ Dx *(dx1+2.*dx2+2.*dx3+dx4)/6.;
                  break;

            case OZAKI :
                  J=cd_m->J_Drift_Function(X);
                  f=cd_m->Drift_Function(X);
                  J1=i(J);
                  R = J1*expm(J*Dx);
                  R=R-J1;
                  U = X + R*f+G*w;
                  break;

            default :
                  break;
            }
      return U;
}

dcovector CD_Simulator::Draw_Init(void)
{
      Discrete_Observed_Model *m=dynamic_cast<Discrete_Observed_Model *>(model);
      dcovector U;
      dcovector X0;
      dsymatrix R0;
      m->Discrete_Observed_Model::Get_Init_Parameters(X0,R0);
      multivariate_normal_draw(U,r,R0);
      return U+X0;
}


dcovector CD_Simulator::Draw_Transition(const dcovector & Xkm1)
{
      int k;
      dcovector Xk=Xkm1;
      for (k=0; k<_a; k++)
                  Xk =draw_state(Xk);
      return Xk;
}

dcovector CD_Simulator::Draw_Observation(const dcovector & Xk)
{
      Discrete_Observed_Model *m=dynamic_cast<Discrete_Observed_Model *>(model);
      dcovector V;
      multivariate_normal_draw(V,r,m->Qv);
      return m->Observation_Function(Xk)+V;
}
long double CD_Simulator::Observation_Density(const dcovector & Yk, const dcovector & Xk)
{
      Discrete_Observed_Model *m=dynamic_cast<Discrete_Observed_Model *>(model);
     return  multivariate_normal_evaluate(Yk-(m->Observation_Function(Xk)),m->Qv);
}

void  CD_Simulator::_update(void)
{
      int N=X.size();
      int k;
      if (N ==0)
            {
                  X.push_back(Draw_Init());
                  Y.push_back(Draw_Observation(X[0]));
            }
      else
            {
                  for (k=0; k<_a; k++)
                        {
                              X.push_back(draw_state(X[X.size()-1]));
                        }
                  Y.push_back(Draw_Observation(X[X.size()-1]));
                        
            }
}


int CD_Simulator::Simulate(const double &T)
{
      dcovector V;
      int Nx = (int)(T/Dx);
      int Ny = (int)(T/Dy);
      int k;
      X.resize(Nx+1);
      Y.resize(Ny+1);

      model->Clear();
      X[0] = Draw_Init();

      Y[0] =Draw_Observation(X[0]);

      for (k=1; k<=Nx; k++)
            {
                  X[k] =draw_state(X[k-1]);
                  if(!(k%_a))
                        {
                              Y[k/_a] =Draw_Observation(X[k]);
                              model->Update();
                        }
            }

      
      return 0;
}

int CD_Simulator::Save_Y(const char *filename)
{
      ofstream file(filename);
      int i,j;
      int N=Y.size();
      double t=0.;

      if(file){
	    for (i=0;i<N;i++)
                  {
                        file << t<<" ";
                        for (j=0; j<Y[0].l; j++)
                              file<<Y[i](j)<<" ";
                        file<<endl;
                        t+=Dy;
	    }
	    file.close();
      }
      
      else
	    {
		  cout<<"ERROR : Pb to access file"<<endl;
		  return 1;
	    }
      return 0;

}

int CD_Simulator::Save_X(const char *filename)
{
      ofstream file(filename);
      int i,j;
      int N=X.size();
      double t=0.;
      
      if(file){
	    for (i=0;i<N;i++)
                  {
                        file << t<<" ";
                        for (j=0; j<X[0].l; j++)
                              file<<X[i](j)<<" ";
                        file<<endl;
                        t+=Dx;
	    }
	    file.close();
      }
      
      else
	    {
		  cout<<"ERROR : Pb to access file"<<endl;
		  return 1;
	    }
      return 0;

}

CD_Simulator_WT::CD_Simulator_WT(void)
{
}
CD_Simulator_WT::CD_Simulator_WT(Continuous_Discrete_Model *cd_m, const int &s,const int &apha, const double &tb, const double &t)
      :CD_Simulator::CD_Simulator(cd_m, s, apha)
{
      TB=tb;
      T=t;
}
      
dcovector CD_Simulator_WT::Draw_Init(void)
{
      int Nx = (int)(T/Dx);
      int Nb = (int)(TB/Dx);
      int k;
      int u;
      if (Xt.size()==0)
            {
                  
                  Xt.resize(Nx+1);

                  Xt[0]=CD_Simulator::Draw_Init();

                  for (k=1; k<=Nb; k++)
                        {
                              Xt[0] =draw_state(Xt[0]);
                        }
 
                  for (k=1; k<=Nx; k++)
                        {
                              Xt[k] =draw_state(Xt[k-1]);
                        }

            }
      u = (int) (gsl_rng_uniform (r)*Nx);
      return Xt[u];
}




LTI_CD_Simulator::LTI_CD_Simulator(void)
{
}
LTI_CD_Simulator::LTI_CD_Simulator(Linear_CD_Model *cd_m, const int &alpha)
{
      _a= alpha;
      Dy=cd_m->Ts;
      Dx =  Dy / _a; 
      model = cd_m;
}
dcovector LTI_CD_Simulator::draw_state(const dcovector &Xkm1)
{
      Linear_CD_Model * m=dynamic_cast<Linear_CD_Model *>(model);
      int N=Xkm1.l;
      int i,j;
      dgematrix P(N,N);
      dcovector W;
      
      dgematrix I(N,N);
      I.identity();

      dgematrix F = expm(m->A * Dx);
      dcovector f = (F - I) * CPPL::i(m->A) *m->B;;
      P.zero();

      // matrix fraction decomposition
      
      dgematrix Ha(2*N,2*N);
      dgematrix CQC = m->C * (m->Qw * t(m->C));
      dgematrix U(2*N,N);
      dgematrix D=P;
      dgematrix E(N,N);
      dgematrix At=t(m->A);
      E.identity();

      for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                  U(i,j) = D(i,j);

      for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                  U(i+N,j) = E(i,j);


      Ha.zero();
      for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                  Ha(i,j) = m->A(i,j);
      for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                  Ha(i,j+N) = CQC(i,j);

      for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                  Ha(i+N,j+N) = -At(i,j);
      U = expm(Ha * Dx) * U;

      for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                  D(i,j) = U(i,j) ;

      for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                  E(i,j) = U(i+N,j);
      P= D*CPPL::i(E);

      if(b)
            W=b(this,r);
      else
            multivariate_normal_draw(W,r,P);
      return F * Xkm1 + f+  W;

}

LTI_CD_Simulator_WT::LTI_CD_Simulator_WT(void)
{
}
LTI_CD_Simulator_WT::LTI_CD_Simulator_WT(Linear_CD_Model *cd_m, const int &apha, const double &tb, const double &t)
      :LTI_CD_Simulator::LTI_CD_Simulator(cd_m, apha)
{
      TB=tb;
      T=t;
}
      
dcovector LTI_CD_Simulator_WT::Draw_Init(void)
{
      int Nx = (int)(T/Dx);
      int Nb = (int)(TB/Dx);
      int k;
      int u;
      if (Xt.size()==0)
            {
                  
                  Xt.resize(Nx+1);

                  Xt[0]=LTI_CD_Simulator::Draw_Init();

                  for (k=1; k<=Nb; k++)
                        {
                              Xt[0] =draw_state(Xt[0]);
                        }
 
                  for (k=1; k<=Nx; k++)
                        {
                              Xt[k] =draw_state(Xt[k-1]);
                        }

            }
      u = (int) (gsl_rng_uniform (r)*Nx);
      return Xt[u];
}

Opt_Simulator::Opt_Simulator(void)
{
}
G_Simulator::G_Simulator(void)
{
}
G_Simulator::G_Simulator(Gaussian_Nonlinear_Model *m)
{
      model=m;
}
dcovector G_Simulator::Draw_Init(void)
{
      Gaussian_Nonlinear_Model *m=dynamic_cast<Gaussian_Nonlinear_Model *>(model);
      dcovector U;
      dcovector X0;
      dsymatrix R0;
      m->Get_Init_Parameters(X0,R0);

      multivariate_normal_draw(U,r,R0);
      U+=X0;

      return U;

}
dcovector G_Simulator::Draw_Transition(const dcovector & Xkm1)
{
      Gaussian_Nonlinear_Model *m=dynamic_cast<Gaussian_Nonlinear_Model *>(model);

      dcovector Xk;
      dcovector W;

      if(b)
            W=b(this,r);
      else
            multivariate_normal_draw(W,r,m->Qw);

      Xk=m->State_Function(Xkm1,W);

      return Xk;

}
dcovector G_Simulator::Draw_Observation(const dcovector & Xk)
{
      Gaussian_Nonlinear_Model *m=dynamic_cast<Gaussian_Nonlinear_Model *>(model);
      dcovector V;
      multivariate_normal_draw(V,r,m->Qv);
      return m->Observation_Function(Xk)+V;

}
long double G_Simulator::Observation_Density(const dcovector & Yk, const dcovector & Xk)
{
      Gaussian_Nonlinear_Model *m=dynamic_cast<Gaussian_Nonlinear_Model *>(model);
      return  multivariate_normal_evaluate(Yk-(m->Observation_Function(Xk)),m->Qv);

}
dcovector G_Simulator::Draw_Optimal(const dcovector & Yk, const dcovector & Xkm1)
{
      Gaussian_Nonlinear_Model *m=dynamic_cast<Gaussian_Nonlinear_Model *>(model);

      dcovector W(m->Qw.n);
      W.zero();

      dgematrix G=m->Jw_State_Function(Xkm1,W);

      dgematrix H=m->J_Observation_Function(Xkm1);
      
      dgematrix M = m->Qw * (t(G) * t(H)) * i(H*G*(m->Qw*t(H*G)) + m->Qv) ;

      dgematrix Q = m->Qw - M * H * G * m->Qw;

      dcovector V;

      dcovector mx=m->State_Function(Xkm1,W) + G * M * (Yk - H * m->State_Function(Xkm1,W));

      multivariate_normal_draw(V,r,Q);  

      dcovector Xs= mx + G * V;

      return Xs;

}
long double G_Simulator::Obs_Optimal_Density(const dcovector & Yk, const dcovector &Xkm1)
{
      Gaussian_Nonlinear_Model *m=dynamic_cast<Gaussian_Nonlinear_Model *>(model);

      dgematrix H=m->J_Observation_Function(Xkm1);

      dcovector W(m->Qw.n);
      W.zero();

      dcovector my = H * m->State_Function(Xkm1,W);

      dgematrix G=m->Jw_State_Function(Xkm1,W);

      dgematrix Q = H * G * (m->Qw * t(H * G)) + m->Qv;
      return multivariate_normal_evaluate(Yk-my,Q);

}

G_Simulator_WT::G_Simulator_WT(void)
{
}
G_Simulator_WT::G_Simulator_WT(Gaussian_Nonlinear_Model *m,const int &nb, const int &n)
      :G_Simulator::G_Simulator(m)
{
      NB=nb;
      N=n;

}
dcovector G_Simulator_WT::Draw_Init(void)
{
      int k;
      int u;
      if (Xt.size()==0)
            {
                  
                  Xt.resize(N);

                  Xt[0]=G_Simulator::Draw_Init();

                  for (k=1; k<=NB; k++)
                        {
                              Xt[0] =Draw_Transition(Xt[0]);
                        }
 
                  for (k=1; k<N; k++)
                        {
                              Xt[k] =Draw_Transition(Xt[k-1]);
                        }

            }
      u = (int) (gsl_rng_uniform (r)*N);
      return Xt[u];
}
