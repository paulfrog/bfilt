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

#include <bfilt/gaussian_model.h>

Model::Model(void)
{
      _k=0;
}
Model::~Model(void)
{
}
int Model::Update(void)
{
      _k++;
      return 0;
}
int Model::Clear(void)
{
      _k=0;
      return 0;
}

int Model::Get_Time(void)
{
      return _k;
}

Discrete_Observed_Model::Discrete_Observed_Model(void)
{
}


Discrete_Observed_Model::~Discrete_Observed_Model(void)
{
}




dgematrix Discrete_Observed_Model::J_Observation_Function(const dcovector & X)
{
      DOM_METHOD p=&Discrete_Observed_Model::Observation_Function;
      return numerical_jacobian(this,p,X);
}

void Discrete_Observed_Model::Get_Init_Parameters(dcovector & mean, dsymatrix &Cov)
{
      mean=X0;
      Cov =R0;
}
Gaussian_Nonlinear_Model::Gaussian_Nonlinear_Model(void)
{
}

Gaussian_Nonlinear_Model::~Gaussian_Nonlinear_Model(void)
{
}


void Gaussian_Nonlinear_Model::Init(void)
{
}


dgematrix Gaussian_Nonlinear_Model::Jx_State_Function(const dcovector & X,const dcovector&  W)
{
      GNM_METHOD_2P p=&Gaussian_Nonlinear_Model::State_Function;
      return numerical_jacobian_p1(this,p,X,W);
}

dgematrix Gaussian_Nonlinear_Model::Jw_State_Function(const dcovector & X,const dcovector&  W)
{
      GNM_METHOD_2P p=&Gaussian_Nonlinear_Model::State_Function;
      return numerical_jacobian_p2(this,p,X,W);
}




void Gaussian_Nonlinear_Model::Get_Linear_Parameters(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &G, dcovector &Xp)
{
      Xp = State_Function(X,W);
      G  = Jw_State_Function(X,W);
      F  = Jx_State_Function(X,W);
}


Gaussian_Linear_Model::Gaussian_Linear_Model(void)
{

}

dcovector Gaussian_Linear_Model::State_Function(const dcovector & X,const dcovector&  W)
{
      return (F*X + f + G * W);
}


dgematrix Gaussian_Linear_Model::Jx_State_Function(const dcovector & X,const dcovector&  W)
{
      return F;
}

dgematrix Gaussian_Linear_Model::Jw_State_Function(const dcovector & X,const dcovector&  W)
{
      return G;
}

dcovector Gaussian_Linear_Model::Get_Mean_Prediction(const dcovector & M)
{
      return F*M + f;
}

dgematrix Gaussian_Linear_Model::Get_Cov_Prediction(const dgematrix & P)
{
      return F * P * t(P) + G * Qw * t(G);
}



dcovector Gaussian_Linear_Model::Observation_Function(const dcovector& X)
{
      return H*X + h;
}

dgematrix Gaussian_Linear_Model::J_Observation_Function(const dcovector & X)
{
      return H;
}


dcovector copy(const dcovector & X, const int & n){
      int i,N=n*X.l;
      dcovector Y(N);
      for (i=0;i<N;i++)
	    Y(i)=X(i%X.l);
      return Y;
}

dsymatrix copy(const dsymatrix & X, const int & a){
      int i,j,k,N=a*X.n;
      dsymatrix Y(N);
      Y.zero();
      for (k=0;k<a;k++)
	    for (j=k*X.n;j<X.n*(k+1);j++){
		  for (i=k*X.n;i<X.n*(k+1);i++)
			Y(i,j)=X(i%X.n,j%X.n);
	    }
      return Y;
}

dgematrix copy(const dgematrix & X, const int & a){
      int i,j,k,N=a*X.n,M=a*X.m;
      dgematrix Y(M,N);
      Y.zero();
      for (k=0;k<a;k++)
	    for (j=k*X.n;j<X.n*(k+1);j++){
		  for (i=k*X.m;i<X.m*(k+1);i++)
			Y(i,j)=X(i%X.m,j%X.n);
	    }
      return Y;
}

int cut(vector<dcovector> & U,const dcovector & X, const int & n){
      int i,j,N=X.l;
      int k=N/n;
      U.resize(n);
      if(!(N%n)){
	    for(j=0;j<n;j++){
		  U[j].resize(k);
		  for(i=0;i<k;i++)
			U[j](i)=X(i+j*k);
	    }    
	    return 0;
      }
      else {
	    cout<< "pb to cut N%n != 0"<<endl;
	    return 1;
      }
}
dcovector paste(const vector<dcovector> & U){
      int j,i,k;
      int n=U.size();
      int l;
      dcovector X;
      if (n!=0){
	    l=U[0].l;
	    X.resize(l*n);
	    k=0;
	    for(i=0;i<n;i++){
		  for(j=0;j<l;j++){
			X(k)=U[i](j);
			k++;
		  }
	    }
      }
      return X;
}

dgematrix paste(const vector<dgematrix> & G){
      int a=G.size();
      dgematrix Y;
      int i,j,k,N,M;
      if(a>0){
	    N=G[0].n;
	    M=G[0].m;
	    Y.resize(a*M,a*N);
	    Y.zero();
	    for (k=0;k<a;k++)
		  for (j=k*M;j<M*(k+1);j++){
			for (i=k*N;i<N*(k+1);i++){
			      Y(j,i)=G[k](j%M,i%N);

			}
		  }
      }
      return Y;
}


Continuous_Discrete_Model::Continuous_Discrete_Model(void)
{
}
Continuous_Discrete_Model::~Continuous_Discrete_Model(void)
{
}

dgematrix Continuous_Discrete_Model::J_Drift_Function(const dcovector & X){
      f_cd_m p=&Continuous_Discrete_Model::Drift_Function;
      return numerical_jacobian(this,p,X);
}

Linear_CD_Model::Linear_CD_Model(void)
{
}

dcovector Linear_CD_Model::Drift_Function(const dcovector & X)
{
      return (A * X + B);
}
  
dgematrix Linear_CD_Model::J_Drift_Function(const dcovector & X)
{
      return A;
}

dgematrix Linear_CD_Model::Diffusion_Function(void)
{
      return C;
}
  
dcovector Linear_CD_Model::Observation_Function(const dcovector& X)
{
      return H*X + h;
}

dgematrix Linear_CD_Model::J_Observation_Function(const dcovector & X)
{
      return H;
}

dcovector Linear_CD_Model::Get_Mean_Prediction(const dcovector & M)
{
      dgematrix EXPA = expm(A*Ts);
      dgematrix I(M.l,M.l);
      I.identity();

      return EXPA * M + (EXPA - I) * CPPL::i(A) *B;
}
dgematrix Linear_CD_Model::Get_Cov_Prediction(const dgematrix & P)
{
      // matrix fraction decomposition
      int i,j;
      int N=P.m;
      dgematrix Ha(2*N,2*N);
      dgematrix CQC = C * (Qw * t(C));
      dgematrix U(2*N,N);
      dgematrix D=P;
      dgematrix E(N,N);
      dgematrix At=t(A);
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
                  Ha(i,j) = A(i,j);
      for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                  Ha(i,j+N) = CQC(i,j);

      for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                  Ha(i+N,j+N) = -At(i,j);
      U = expm(Ha * Ts) * U;
      
      for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                  D(i,j) = U(i,j) ;

      for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                  E(i,j) = U(i+N,j);

      return D*CPPL::i(E);
}




//  -----------------------------               Discrete_Approximation_CD_Model


Discrete_Approximation_CD_Model::Discrete_Approximation_CD_Model(void): Gaussian_Nonlinear_Model(){
      cd_model=NULL;
      alpha=0;
}

Discrete_Approximation_CD_Model::~Discrete_Approximation_CD_Model(void)
{
}

Discrete_Approximation_CD_Model::Discrete_Approximation_CD_Model(Continuous_Discrete_Model *m)
{
      dcovector x0;
      dsymatrix r0;
      alpha=1;
      cd_model=m;
      cd_model->Get_Init_Parameters(x0,r0);

      X0=x0;
      R0=r0;
      Qw=cd_model->Qw * cd_model->Ts;
      Qv=cd_model->Qv;
}
Discrete_Approximation_CD_Model::Discrete_Approximation_CD_Model(Continuous_Discrete_Model *m,const int &a)
{
      double T;
      dsymatrix Q;
      dcovector x0;
      dsymatrix r0;
      
      alpha=a;
      cd_model=m;
      Q=cd_model->Qw;
      T= cd_model->Ts/((double)alpha);

      cd_model->Get_Init_Parameters(x0,r0);

      X0=x0;
      R0=r0;

      Q*=T;
      Qw=copy(Q,alpha);
      Qv=cd_model->Qv;

}

void Discrete_Approximation_CD_Model::Set_Alpha(const int & a)
{
      double T;
      dsymatrix Q;
  
      alpha=a;
      Q=cd_model->Qw;

      T= cd_model->Ts/((double)alpha);

      Q*=T;

      Qw=copy(Q,alpha);
}

int Discrete_Approximation_CD_Model::Get_Alpha(void)
{
      return alpha;
}

void Discrete_Approximation_CD_Model::Get_Linear_Parameters(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &J, dcovector &Xp){
      int i;
      vector<dcovector> Xk(alpha+1);
      vector<dcovector> Uw;
      vector<dgematrix> Jx(alpha);
      vector<dgematrix> Jw(alpha);
      vector<dgematrix> L(alpha);
      vector<dgematrix> G(alpha);


      cut(Uw,W,alpha);
      Xk[0]=X;
      F.resize(X.l,X.l);
      F.identity();

      for(i=0;i<alpha;i++)
            {
                  Get_Linear_Scheme(Xk[i],Uw[i],Jx[i],Jw[i],Xk[i+1]);
                  F*=Jx[i];
            }
      Xp= Xk[alpha];

      

      L[0].resize(X.l,X.l);
      L[0].identity();

      for(i=1; i<alpha; i++)
            L[i] =  L[i-1] * Jx[alpha - i];

      for(i=0; i<alpha; i++)
            G[i] =  L[alpha-i-1] * Jw[i];

      J=cat(G); 
}


void Discrete_Approximation_CD_Model::Init(void){
      cd_model->Init();
      dcovector x0;
      dsymatrix r0;

      cd_model->Get_Init_Parameters(x0,r0);
      X0=x0;
      R0=r0;
}
dcovector Discrete_Approximation_CD_Model::State_Function(const dcovector & X,const dcovector&  W)
{
      int i;
      dcovector Xk;
      vector<dcovector> Uw;

      cut(Uw,W,alpha);

      Xk=X;
      for(i=0;i<alpha;i++)
	    Xk=Scheme(Xk,Uw[i]);

      return Xk;
}

dgematrix Discrete_Approximation_CD_Model::Jx_State_Function(const dcovector & X,const dcovector&  W)
{
      int i;
      vector<dcovector> Xk(alpha);
      vector<dcovector> Uw;
      dgematrix J(X.l,X.l);
      cut(Uw,W,alpha);
      
      Xk[0]=X;
      for(i=1;i<alpha;i++)
            Xk[i]=Scheme(Xk[i-1],Uw[i-1]);
      
      J.identity();
      
      for(i=0; i<alpha; i++)
            J *= Jx_Scheme(Xk[i],Uw[i]); 

      return J;
}

dgematrix Discrete_Approximation_CD_Model::Jw_State_Function(const dcovector & X,const dcovector&  W)
{
      dgematrix J(X.l,W.l);
      int i;
      vector<dcovector> Xk(alpha);
      vector<dcovector> Uw;
      vector<dgematrix> Jx(alpha);
      vector<dgematrix> Jw(alpha);
      vector<dgematrix> L(alpha);
      vector<dgematrix> G(alpha);


      cut(Uw,W,alpha);
      
      Xk[0]=X;
      for(i=1;i<alpha;i++)
            Xk[i]=Scheme(Xk[i-1],Uw[i-1]);
      
      for(i=0; i<alpha; i++)
            Jx[i] = Jx_Scheme(Xk[i],Uw[i]); 
      
      for(i=0; i<alpha; i++)
            Jw[i] = Jw_Scheme(Xk[i],Uw[i]); 

      L[0].resize(X.l,X.l);
      L[0].identity();

      for(i=1; i<alpha; i++)
            L[i] =  L[i-1] * Jx[alpha - i];

      for(i=0; i<alpha; i++)
            G[i] =  L[alpha-i-1] * Jw[i];

      J=cat(G); 
      
      return J;

            
}

  
dcovector Discrete_Approximation_CD_Model::Observation_Function(const dcovector &X)
{
      return cd_model->Observation_Function(X);
}

dgematrix Discrete_Approximation_CD_Model::J_Observation_Function(const dcovector& X)
{
      return cd_model->J_Observation_Function(X);
}



//   ------------------------          EULER  

Euler_CD_Model::Euler_CD_Model(void)
{

}
Euler_CD_Model::Euler_CD_Model(Continuous_Discrete_Model *m, const int & a):
      Discrete_Approximation_CD_Model(m,a)
{
}

dcovector Euler_CD_Model::Scheme(const dcovector & X,const dcovector & W)
{
      double delta = cd_model->Ts /( (double) alpha );
      return (X + delta*cd_model->Drift_Function(X) + cd_model->Diffusion_Function()*W);
}

dgematrix Euler_CD_Model::Jx_Scheme(const dcovector &X, const dcovector &W)
{
      dgematrix L(X.l);

      L.identity();

      return 	(L + cd_model->Ts/( (double) alpha ) * cd_model->J_Drift_Function(X));

}

dgematrix Euler_CD_Model::Jw_Scheme(const dcovector &X, const dcovector &W)
{
      return cd_model->Diffusion_Function();
}

void Euler_CD_Model::Get_Linear_Scheme(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &J, dcovector &Xp)
{
      double delta = cd_model->Ts /( (double) alpha );

      dgematrix L(X.l);

      L.identity();

      Xp = (X + delta*cd_model->Drift_Function(X) + cd_model->Diffusion_Function()*W);
      
      F = (L + cd_model->Ts/( (double) alpha ) * cd_model->J_Drift_Function(X));
      
      J = cd_model->Diffusion_Function();
}

//       --------------------        SRK4  :

SRK4_CD_Model::SRK4_CD_Model(void)
{
}
SRK4_CD_Model::SRK4_CD_Model(Continuous_Discrete_Model *m,const int & a):
      Discrete_Approximation_CD_Model(m,a)
{
}


dcovector SRK4_CD_Model::Scheme(const dcovector & X,const dcovector & W){

      dcovector U(X.l), dx1(X.l), dx2(X.l), dx3(X.l), dx4(X.l);

      double delta = cd_model->Ts /( (double) alpha );  

      dgematrix G=cd_model->Diffusion_Function();

      dx1 =  cd_model->Drift_Function(X) + G*W/delta;

      dx2 =  cd_model->Drift_Function(X + 0.5*delta *dx1)+G*W/delta;

      dx3 =  cd_model->Drift_Function(X + 0.5*delta *dx2)+G*W/delta;

      dx4 =  cd_model->Drift_Function(X + delta *dx3) + G*W/delta; 
  
      U=X+ delta *(dx1+2.*dx2+2.*dx3+dx4)/6.;
  
      return U;
}

dgematrix SRK4_CD_Model::Jx_Scheme(const dcovector &X, const dcovector &W)
{

      dgematrix J1,J2,J3,J4,J,I(X.l,X.l);
      I.identity();
      dcovector U(X.l), dx1(X.l), dx2(X.l), dx3(X.l), dx4(X.l);

      double delta = cd_model->Ts/( (double) alpha );  

      dgematrix G=cd_model->Diffusion_Function();

      dx1 = cd_model->Drift_Function(X) + G*W/delta;

      dx2 = cd_model->Drift_Function(X + 0.5*delta*dx1) + G*W/delta;

      dx3 = cd_model->Drift_Function(X + 0.5*delta*dx2) + G*W/delta;

      dx4 = cd_model->Drift_Function(X +     delta*dx3) + G*W/delta; 
  
      J1= cd_model->J_Drift_Function(X);
      J2= cd_model->J_Drift_Function(X + 0.5*delta*dx1) * (I+0.5*delta*J1);
      J3= cd_model->J_Drift_Function(X + 0.5*delta*dx2) * (I+0.5*delta*J2);
      J4= cd_model->J_Drift_Function(X +     delta*dx3) * (I+    delta*J3);


      J=I+ delta * (J1+2.*J2+2.*J3+J4)/6.;
  
      return J;
}

dgematrix SRK4_CD_Model::Jw_Scheme(const dcovector &X, const dcovector &W)
{
      dgematrix J1,J2,J3,J4,J;
      
      dcovector U(X.l), dx1(X.l), dx2(X.l), dx3(X.l), dx4(X.l);

      double delta = cd_model->Ts/( (double) alpha );  

      dgematrix G=cd_model->Diffusion_Function();


      dx1 = cd_model->Drift_Function(X)+ G*W/delta;;

      dx2 = cd_model->Drift_Function(X + 0.5*delta*dx1) + G*W/delta;;

      dx3 = cd_model->Drift_Function(X + 0.5*delta*dx2) + G*W/delta;;

      dx4 = cd_model->Drift_Function(X +     delta*dx3) + G*W/delta;;

      J1= G/delta;
      J2= cd_model->J_Drift_Function(X + 0.5*delta*dx1) * (0.5*delta*J1)  + G/delta;
      J3= cd_model->J_Drift_Function(X + 0.5*delta*dx2) * (0.5*delta*J2)  + G/delta;
      J4= cd_model->J_Drift_Function(X +     delta*dx3) * (    delta*J3)  + G/delta;


      J=(J1+2.*J2+2.*J3+J4)*delta/6.;

      return J;
}

void SRK4_CD_Model::Get_Linear_Scheme(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &Jw, dcovector &Xp)
{
      double delta = cd_model->Ts /( (double) alpha );
      dcovector U(X.l), dx1(X.l), dx2(X.l), dx3(X.l), dx4(X.l);
      dgematrix J1,J2,J3,J4,I(X.l,X.l);
      I.identity();
      dgematrix G=cd_model->Diffusion_Function();

      // compute the prediction mean Xp
      dx1 =  cd_model->Drift_Function(X) + G*W/delta;
      
      dx2 =  cd_model->Drift_Function(X + 0.5*delta *dx1)+G*W/delta;

      dx3 =  cd_model->Drift_Function(X + 0.5*delta *dx2)+G*W/delta;

      dx4 =  cd_model->Drift_Function(X + delta *dx3) + G*W/delta; 
  
      Xp  =  X+ delta *(dx1+2.*dx2+2.*dx3+dx4)/6.;

      // compute Jx 

      J1= cd_model->J_Drift_Function(X);

      J2= cd_model->J_Drift_Function(X + 0.5*delta*dx1) * (I+0.5*delta*J1);
     
      J3= cd_model->J_Drift_Function(X + 0.5*delta*dx2) * (I+0.5*delta*J2);
     
      J4= cd_model->J_Drift_Function(X +     delta*dx3) * (I+    delta*J3);

           
      F=I+ delta * (J1+2.*J2+2.*J3+J4)/6.;
     

      // compude Jw
      J1= G/delta;
     
      J2= cd_model->J_Drift_Function(X + 0.5*delta*dx1) * (0.5*delta*J1)  + G/delta;
     
      J3= cd_model->J_Drift_Function(X + 0.5*delta*dx2) * (0.5*delta*J2)  + G/delta;
     
      J4= cd_model->J_Drift_Function(X +     delta*dx3) * (    delta*J3)  + G/delta;
     
     
      Jw=(J1+2.*J2+2.*J3+J4)*delta/6.;

}
//       --------------------        Heun  :

Heun_CD_Model::Heun_CD_Model(void)
{
}
Heun_CD_Model::Heun_CD_Model(Continuous_Discrete_Model *m,const int & a):
      Discrete_Approximation_CD_Model(m,a)
{
}


dcovector Heun_CD_Model::Scheme(const dcovector & X,const dcovector & W){

      dcovector U(X.l), dx1(X.l), dx2(X.l);

      double delta = cd_model->Ts /( (double) alpha );  

      dgematrix G=cd_model->Diffusion_Function();

      dx1 = cd_model->Drift_Function(X) + G*W/delta;

      dx2 = cd_model->Drift_Function(X + dx1*delta) + G*W/delta; 
  
      U=X+ (dx1+dx2)*delta/2.;
  
      return U;
}

dgematrix Heun_CD_Model::Jx_Scheme(const dcovector &X, const dcovector &W)
{

      dgematrix J1,J2,J,I(X.l,X.l);
      I.identity();
      dcovector U(X.l), dx1(X.l), dx2(X.l);

      double delta = cd_model->Ts/( (double) alpha );  

      dgematrix G=cd_model->Diffusion_Function();

      dx1 = cd_model->Drift_Function(X) + G*W/delta;

      dx2 = cd_model->Drift_Function(X +     delta*dx1) + G*W/delta; 
  
      J1= cd_model->J_Drift_Function(X);
      J2= cd_model->J_Drift_Function(X +     delta*dx1) * (I+    delta*J1);


      J=I+ delta * (J1+J2)/2.;
  
      return J;
}

dgematrix Heun_CD_Model::Jw_Scheme(const dcovector &X, const dcovector &W)
{
      dgematrix J1,J2,J;
      
      dcovector U(X.l), dx1(X.l), dx2(X.l);

      double delta = cd_model->Ts/( (double) alpha );  

      dgematrix G=cd_model->Diffusion_Function();


      dx1 = cd_model->Drift_Function(X) + G*W/delta;

      dx2 = cd_model->Drift_Function(X + delta*dx1)  + G*W/delta;

      J1= G/delta;

      J2= cd_model->J_Drift_Function(X +     delta*dx1) * (    delta*J1)  + G/delta;


      J=(J1+J2)*delta/2.;

      return J;
}

void Heun_CD_Model::Get_Linear_Scheme(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &J, dcovector &Xp)
{
      dgematrix J1,J2,I(X.l,X.l);
      I.identity();
      dcovector U(X.l), dx1(X.l), dx2(X.l);

      double delta = cd_model->Ts/( (double) alpha );  

      dgematrix G=cd_model->Diffusion_Function();

      dx1 = cd_model->Drift_Function(X) + G*W/delta;

      dx2 = cd_model->Drift_Function(X + delta*dx1)  + G*W/delta;

      Xp=X+ (dx1+dx2)*delta/2.;

      J1= cd_model->J_Drift_Function(X);

      J2= cd_model->J_Drift_Function(X +     delta*dx1) * (I+    delta*J1);


      F=I+ delta * (J1+J2)/2.;

      J1= G/delta;

      J2= cd_model->J_Drift_Function(X +     delta*dx1) * (    delta*J1)  + G/delta;


      J=(J1+J2)*delta/2.;
}

//       --------------------        OZAKI  :

Ozaki_CD_Model::Ozaki_CD_Model(void)
{
}
Ozaki_CD_Model::Ozaki_CD_Model(Continuous_Discrete_Model *m,const int & a):
      Discrete_Approximation_CD_Model(m,a)
{
}


dcovector Ozaki_CD_Model::Scheme(const dcovector & X,const dcovector & W){
      double delta = cd_model->Ts /( (double) alpha );  

      dgematrix G=cd_model->Diffusion_Function();
      dgematrix J=cd_model->J_Drift_Function(X);
      dcovector f=cd_model->Drift_Function(X);

      dgematrix J1=i(J);
      dgematrix R = J1*expm(J*delta);
      R=R-J1;
	
      dcovector U;
  
      return U = X + R*f+G*W;
}

dgematrix Ozaki_CD_Model::Jx_Scheme(const dcovector &X, const dcovector &W)
{

      double delta = cd_model->Ts /( (double) alpha );  

      dgematrix J=cd_model->J_Drift_Function(X);

      dgematrix R = expm(J*delta);
  
      return R;

}

dgematrix Ozaki_CD_Model::Jw_Scheme(const dcovector &X, const dcovector &W)
{
      return cd_model->Diffusion_Function();
}

void Ozaki_CD_Model::Get_Linear_Scheme(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &G, dcovector &Xp)
{
      double delta = cd_model->Ts /( (double) alpha );  
     
      G=cd_model->Diffusion_Function();
      dgematrix J=cd_model->J_Drift_Function(X);
      dcovector f=cd_model->Drift_Function(X);
      
      dgematrix J1=i(J);
      F = expm(J*delta);
      dgematrix R = J1*F;
      R=R-J1;
      
      Xp = X + R*f+G*W;
 
}

dgematrix cat(const vector<dgematrix> &V)
{
      int S= V.size();
      int L = V[0].m;
      int C = V[0].n;
      int i,j,k;
      dgematrix J(L,C*S);

      for (i=0; i<L; i++)
            for (j=0; j<C*S; j++)
                  J(i,j) = V[(int)(j/C)](i%L,j%C);
           
      return J;
      
}
