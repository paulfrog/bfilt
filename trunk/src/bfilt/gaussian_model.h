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

///
/// @file   gaussian_model.h
/// @author paul <paul.frogerais@univ-rennes1.fr>
/// @date   Fri Sep 12 18:35:22 2008
/// 
/// @brief  Implementation of gaussian non linear models
/// 
/// 
///

#ifndef GAUSS__NLIN__MODEL
#define GAUSS__NLIN__MODEL



#include <bfilt/filter_tools.h>

/// @brief The class of time varying-models
/// 
///
class Model
{
protected :
      /// @brief The time
      int _k;
public :
      Model(void);
      /// @brief The Destructor
      ///
      ///
      /// @return 
      ///
      virtual ~Model(void);

      /// @brief Update the time
      ///
      ///
      /// @return 0 if it's Ok;
      ///
      int Update(void);

      /// @brief Set the time to 0
      ///
      ///
      /// @return 0 if it's Ok;
      ///
      int Clear(void);

      /// @brief Get The current time
      ///
      ///
      /// @return 
      ///      
      int Get_Time(void);
};

/// @brief Class of discretely observed model
/// 
/// The output \f$Y_k\f$ is a discrete form of the hidden state.
/// The init state is gaussian \f$ \sim \mathcal{N}(X0,R0)\f$.
/// The state  and observation noises \f$ Wk, V_k \f$ are zero-mean gaussians processes.
/// Their respective covariances are  \f$ Q_w \f$ and \f$ Q_v\f$  .

class Discrete_Observed_Model : public Model
{
protected :
      /// @brief The covariance matrix of p(X0)
      dsymatrix R0;

      /// @brief The mean of p(X0)
      dcovector X0;
public :
      /// @brief The covariance matrix of state noise 
      dsymatrix Qw;

      /// @brief The covariance matrix of observation noise
      dsymatrix Qv;


      /// @brief The  Constructor
      ///
      ///
      Discrete_Observed_Model(void);

      /// @brief The Destructor
      ///
      ///
      /// @return 
      ///
      virtual ~Discrete_Observed_Model(void);

      /// @brief The observation Yk=H(Xk) + Vk
      ///
      /// @param X The state at k
      ///
      /// @return The observation at k
      ///
      virtual dcovector Observation_Function(const dcovector& X)=0;

      /// @brief the jacobian of the observation function
      /// evaluate at X
      /// @param X 
      ///
      /// @return The jacobian matrix
      ///
      virtual dgematrix J_Observation_Function(const dcovector & X);

      /// @brief Return the first an second moment of the initial law p(X0)
      ///
      /// @param mean The mean X0
      /// @param Cov  The Covariance R0
      ///
      virtual void Get_Init_Parameters(dcovector & mean, dsymatrix &Cov);
      

};

/// @brief Gaussian Nonlinear Model
/// The state : X(k) = F (Xk-1, Wk)
/// The Observation Y(k) = H (X(k)) + V
///

class Gaussian_Nonlinear_Model: public Discrete_Observed_Model
{
public :

      


      Gaussian_Nonlinear_Model(void);
      virtual ~Gaussian_Nonlinear_Model(void);

      /// @brief Init The model if needed
      ///
      /// 
      virtual void  Init(void); 

      /// @brief The state Xk=F(Xk-1,Wk)
      ///
      /// @param X The state at k-1
      /// @param W The Noise 
      ///
      /// @return The state at k
      ///
      virtual dcovector State_Function(const dcovector & X,const dcovector&  W)=0;

      /// \brief the X jacobian of the State function
      /// evaluate at X,W
      /// @param X evaluate at X
      /// @param W evalate at W
      ///
      /// @return The jacobian matrix
      ///
      virtual dgematrix Jx_State_Function(const dcovector & X,const dcovector&  W);

      /// \brief the W jacobian of the State function
      /// evaluate at X,W
      /// @param X evaluate at X
      /// @param W evalate at W
      ///
      /// @return The jacobian matrix
      ///
      virtual dgematrix Jw_State_Function(const dcovector & X,const dcovector&  W);

      /// \brief computed linearized parameter for EKF in X,W
      ///
      /// @param X The state value
      /// @param W The noise value
      /// @param F The jacobian of f(X,W) in X
      /// @param G The jacobian in f(X,W) in W
      /// @param Xp The prediction Xp = f(X,W)
      ///
      virtual void Get_Linear_Parameters(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &G, dcovector &Xp);
};

/// @brief Gaussian Linear Model : 
///
/// The state : X(k) = F X(k-1) + f + G * Wk
/// The Observation Y(k) = H X(k) + h + V
///
///

class Gaussian_Linear_Model: public Gaussian_Nonlinear_Model
{
public :
      dgematrix F;
      dgematrix G;
      dcovector f;
      dcovector h;
      dgematrix H;

public :
      Gaussian_Linear_Model(void);
      dcovector State_Function(const dcovector & X,const dcovector&  W);
      dgematrix Jx_State_Function(const dcovector & X,const dcovector&  W);
      dgematrix Jw_State_Function(const dcovector & X,const dcovector&  W);
      dcovector Get_Mean_Prediction(const dcovector & M);
      dgematrix Get_Cov_Prediction(const dgematrix & P);      
      dcovector Observation_Function(const dcovector& X);
      dgematrix J_Observation_Function(const dcovector & X);


};


/// \brief Continuous Discrete Model:
/// The continuous state : \f$ dX(t) = F(X)dt + G()*d\beta \f$
/// The discrete Observation \f$ Yk = H (X(tk)) + Vk \f$

class Continuous_Discrete_Model : public Discrete_Observed_Model
{
public :

      /// \brief The sampling periode Ts=tk - tk-1
      double Ts;

      /// The constructor
      ///
      ///
      Continuous_Discrete_Model(void);

      /// The destructor
      ///
      ///
      virtual ~Continuous_Discrete_Model(void);

      /// \brief the drift function of dX(t) = F(X)dt + G(X)*dW
      ///
      /// @param X the state. 
      ///
      /// @return f(X)
      ///
      virtual dcovector Drift_Function(const dcovector & X)=0;
  
      /// \brief the jacobian of the drift function
      /// evaluate at X
      /// @param X 
      ///
      /// @return the jacobian matrix
      ///
      virtual dgematrix J_Drift_Function(const dcovector & X);

      /// \brief the diffusion function
      ///
      /// @param X the state 
      ///
      /// @return G(X). 
      ///      
      virtual dgematrix Diffusion_Function(void)=0;
    
      /// Initialized CD model
      ///
      ///
      virtual void Init(void){};

};
typedef  dcovector (Continuous_Discrete_Model::*f_cd_m)(const dcovector &  x);


///
/// @brief Linear continuous discrete model class
/// of the form dx = AX dt + Bdt + CdW
/// Y_k = HX(t_k) + h + V_k
///

class Linear_CD_Model: public Continuous_Discrete_Model
{
public :
      dgematrix A;
      dcovector B;
      dgematrix C;
      dcovector h;
      dgematrix H;

      Linear_CD_Model(void);

      dcovector Drift_Function(const dcovector & X);
  
      dgematrix J_Drift_Function(const dcovector & X);

      dgematrix Diffusion_Function(void);

      dcovector Observation_Function(const dcovector& X);

      dgematrix J_Observation_Function(const dcovector & X);

      dcovector Get_Mean_Prediction(const dcovector & M);

      dgematrix Get_Cov_Prediction(const dgematrix & P);      

      virtual void Init(void){};

};

/// @brief the continuous state equation is discretly approximate
/// by X(tk) = f'(X(tk-1),Wk) 
///
///
class Discrete_Approximation_CD_Model: public  Gaussian_Nonlinear_Model
{
protected :

      /// \brief the continuous discrete model
      Continuous_Discrete_Model *cd_model;

      /// \brief  the resolution of the discrete step
      /// Td = Ts * a (Ts = sample duration of discrete observation)
      int alpha;
	
public:


      Discrete_Approximation_CD_Model(void);
      virtual ~Discrete_Approximation_CD_Model(void);
      Discrete_Approximation_CD_Model(Continuous_Discrete_Model *m);

      ///  \brief the constructor
      ///
      /// @param m the CD model
      /// @param a the resolution of the discrete step Td = Ts * a (Ts = sample duration of discrete observation)
      ///
      ///
      Discrete_Approximation_CD_Model(Continuous_Discrete_Model *m, const int & a);

      dcovector State_Function(const dcovector & X,const dcovector&  W);
      dgematrix Jx_State_Function(const dcovector & X,const dcovector&  W);

      void Get_Linear_Parameters(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &G, dcovector &Xp);

      /// \brief Get the Linearized parameters Scheme in X,W
      ///
      /// @param X The state value
      /// @param W The noise value
      /// @param F The jacobian of f(X,W) in X
      /// @param G The jacobian in f(X,W) in W
      /// @param Xp The prediction Xp = f(X,W)
      ///
      virtual  void Get_Linear_Scheme(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &J, dcovector &Xp)=0;

      dgematrix Jw_State_Function(const dcovector & X,const dcovector&  W);

      dcovector Observation_Function(const dcovector& X);

      virtual dgematrix J_Observation_Function(const dcovector& X);

      virtual  dcovector Scheme(const dcovector &X, const dcovector &W)=0; 

      virtual  dgematrix Jx_Scheme(const dcovector &X, const dcovector &W)=0; 

      virtual  dgematrix Jw_Scheme(const dcovector &X, const dcovector &W)=0; 

      virtual void Init(void);

      void Set_Alpha(const int & a);

      int Get_Alpha(void);
};


/// \brief continuous discret model: 
/// the state SDE is discretly approximate by an Euler method 
///
class Euler_CD_Model : public Discrete_Approximation_CD_Model
{

public :
      Euler_CD_Model(void);

      Euler_CD_Model(Continuous_Discrete_Model *m, const int & a);

      dcovector Scheme(const dcovector & X,const dcovector & W);

      void Get_Linear_Scheme(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &J, dcovector &Xp);

      dgematrix Jx_Scheme(const dcovector &X, const dcovector &W);

      dgematrix Jw_Scheme(const dcovector &X, const dcovector &W);
};




/// \brief continuous discret model: 
/// the state SDE is discretly approximate by an Sstochastic runge kutta  method 
///

class SRK4_CD_Model :  public Discrete_Approximation_CD_Model
{

public :
      SRK4_CD_Model(void);

      SRK4_CD_Model(Continuous_Discrete_Model *m,const int & a);

      dcovector Scheme(const dcovector & X,const dcovector & W);

      void Get_Linear_Scheme(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &J, dcovector &Xp);

      dgematrix Jx_Scheme(const dcovector &X, const dcovector &W);

      dgematrix Jw_Scheme(const dcovector &X, const dcovector &W);
};

/// \brief continuous discret model: 
/// the state SDE is discretly approximate by an Sstochastic Heun  method 
///
class Heun_CD_Model : public Discrete_Approximation_CD_Model
{

public :
      Heun_CD_Model(void);

      Heun_CD_Model(Continuous_Discrete_Model *m,const int & a);

      dcovector Scheme(const dcovector & X,const dcovector & W);

      dgematrix Jx_Scheme(const dcovector &X, const dcovector &W);

      dgematrix Jw_Scheme(const dcovector &X, const dcovector &W);

      void Get_Linear_Scheme(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &J, dcovector &Xp);
};

/// \brief continuous discret model: 
/// the state SDE is discretly approximate by Ozaki  method 
///
class Ozaki_CD_Model : public Discrete_Approximation_CD_Model
{

public :
      Ozaki_CD_Model(void);

      Ozaki_CD_Model(Continuous_Discrete_Model *m,const int & a);

      dcovector Scheme(const dcovector & X,const dcovector & W);

      dgematrix Jx_Scheme(const dcovector &X, const dcovector &W);

      dgematrix Jw_Scheme(const dcovector &X, const dcovector &W);
      
      void Get_Linear_Scheme(const dcovector &X,const dcovector &W,dgematrix &F,dgematrix &J, dcovector &Xp);
};





//       ----     useful functions :

dcovector paste(const vector<dcovector> & U);

int cut(vector<dcovector> & U,const dcovector & X, const int & n);

dsymatrix copy(const dsymatrix & X, const int & a);

dcovector copy(const dcovector & X, const int & n);

dgematrix copy(const dgematrix & X, const int & a);

dgematrix paste(const vector<dgematrix> & G);


dgematrix cat(const vector<dgematrix> &V);



typedef  dcovector (Discrete_Observed_Model::*DOM_METHOD)(const dcovector &  x);
typedef  dcovector (Gaussian_Nonlinear_Model::*GNM_METHOD)(const dcovector &  x);
typedef  dcovector (Gaussian_Nonlinear_Model::*GNM_METHOD_2P)(const dcovector &  x,const dcovector &  w);

template<class T,typename F>
dgematrix numerical_jacobian(T * m, F function, const dcovector & X){
      int i,j,N=X.l;
      dcovector Ji,dx(N);
      dgematrix J;
      double eps=1e-10,x;
      dx.zero();
  
      dx(0)=eps;	
      Ji=((m->*function)(X+dx)-(m->*function)(X))/eps;
      J.resize(Ji.l,N);
      for(j=0;j<Ji.l;j++)
	    J(j,0)=Ji(j);
      dx.zero();
      for(i=1;i<N;i++){
	    dx(i)=eps;	
	    Ji=((m->*function)(X+dx)-(m->*function)(X))/eps;
    
	    for(j=0;j<Ji.l;j++)
		  J(j,i)=Ji(j);
    
	    dx.zero();
      }
      return J;	
}


template<class T,typename F>
dgematrix numerical_jacobian_p1(T * m, F function, const dcovector & X,const dcovector &W){
      int i,j,N=X.l;
      dcovector Ji,dx(N);
      dgematrix J;
      double eps=1e-10,x;
      dx.zero();
  
      dx(0)=eps;	
      Ji=((m->*function)(X+dx,W)-(m->*function)(X,W))/eps;
      J.resize(Ji.l,N);
      for(j=0;j<Ji.l;j++)
	    J(j,0)=Ji(j);
      dx.zero();
      for(i=1;i<N;i++){
	    dx(i)=eps;	
	    Ji=((m->*function)(X+dx,W)-(m->*function)(X,W))/eps;
    
	    for(j=0;j<Ji.l;j++)
		  J(j,i)=Ji(j);
    
	    dx.zero();
      }
      return J;	
}


template<class T,typename F>
dgematrix numerical_jacobian_p2(T * m, F function, const dcovector & X,const dcovector &W){
      int i,j,N=W.l;
      dcovector Ji,dw(N);
      dgematrix J;
      double eps=1e-10,w;
      dw.zero();
  
      dw(0)=eps;	
      Ji=((m->*function)(X,W+dw)-(m->*function)(X,W))/eps;
      J.resize(Ji.l,N);
      for(j=0;j<Ji.l;j++)
	    J(j,0)=Ji(j);
      dw.zero();
      for(i=1;i<N;i++){
	    dw(i)=eps;	
	    Ji=((m->*function)(X,W + dw)-(m->*function)(X,W))/eps;
    
	    for(j=0;j<Ji.l;j++)
		  J(j,i)=Ji(j);
    
	    dw.zero();
      }
      return J;	
}



#endif

