/// @mainpage
///
/// @section description Description 
///
/// BFilt is  a multi-platform and open-source  C++ bayesian filtering
/// library.  It contains  useful  and classical  algorithms in  state
/// estimation of  hidden markov models.  So you can  easily construct
/// discrete-disrete (DD) and  continuous-discrete (CD) models (linear
/// or nonlinear)  for filtering (Kalman, EKF,  UKF, particle filters,
/// ...)  and simulation methods.  Indeed, markovian  model simulators
/// can  be used  for  particle  filters. Libraries  such  as BFL  and
/// Bayes++ consider only discrete-discrete filtering. With BFilt, you
/// can easily construct  your own CD or DD  models for filtering. For
/// CD models  stochastic discretization methods  (Euler, Runge Kutta,
/// Local  linearization,  Heun)  are  implemented in  simulation  and
/// filtering.
///
/// @section implementation Dependances
///
/// LAPACK  and  CPPLAPACK  libraries  are  used  for  linear  algebra
/// operations.  For best  performances it  is recommended  to compile
/// yourself  the LAPPACK  libraries  with ATLAS.  The Gnu  Scientific
/// Library  (GSL)  achieves   random  drawing  in  simulators.  These
/// open-source and multi-platform  libraries must be installed before
/// install BFilt.
///
/// @section Installation
/// Go to the bin directory
/// @code
/// cd BFilt/bin
/// @endcode 
/// Run Cmake (>2.6)
/// @code
/// cmake ../src
/// @endcode
/// Compile Bfilt
/// @code 
/// make
/// @endcode 
/// Install BFilt in /usr/local/lib or /usr/local/inlcude 
/// @code 
/// make install
/// @endcode
/// If you want to change the default install directory you can type
/// @code 
/// ccmake 
/// @endcode 
/// and change CMAKE_INSTALL_PREFIX 
/// @section auteur Auteur
///
/// @author paul <paul.frogerais@univ-rennes1.fr>
/// @date   Fri Sep 12 18:34:36 2008
///
///
/// @file   filter.h
/// @author paul <paul.frogerais@univ-rennes1.fr>
/// @date   Fri Sep 12 18:34:36 2008
/// 
/// @brief  Abstract classes of filters
/// 
/// 
///

#ifndef __FILTER__H
#define  __FILTER__H

#include <bfilt/gaussian_model.h>

///  @brief Abstract class of all filters.
/// 
///  Filters calculate recursively an estimation \f$ \hat{X}_{k|k} \f$
///  of the STATE \f$ X_{k} \f$ of a hidden markov model (Model) given
///  observations \f$ Y_{0:k} \f$. \par
///
///  They    compute    also    recursively   the    likelihood    
///  \f$ p_{Y_{0:N}}(y_{0:N}) \f$.

class Filter
{
protected :
      /// The likelihood \f$ p_{Y_{0:N}}(y_{0:N}) \f$
      double Likelihood;
public :

      /// The hidden markov  model 
      Model *model;

      /// @brief \f$ \{ \hat{X}_{k|k} ,k=0,...N \} \f$
      /// 
      /// The estimated trajectory of the state
      vector<dcovector> X;
      
      /// A constructor
      ///
      ///
      Filter(void);

      /// The destructor
      ///
      ///
      virtual ~Filter(void);

      /// Perform an estimation step with a new observation
      ///
      /// @param Y The new observed sample
      ///
      /// @return 0 if everything is ok
      ///
      int Update(const dcovector & Y);

      /// Perform a trajectory state estimation given 
      /// a sequence \f$ y_{0:N}\f$
      ///
      /// @param Y The sequence
      ///
      /// @return 0 if everything is ok
      ///
      int Filtering(const vector<dcovector> & Y);

      /// Evaluate the current estimation of the state
      ///
      ///
      /// @return \f$ \hat{X}_{k|k} \f$
      ///
      virtual dcovector Expected_Get(void)=0;

      /// To init the filter at k=0
      ///
      ///
      int Init(void); 

      /// Return the current likelihood \f$ p_{Y_{0:k}}(y_{0:k}) \f$
      ///
      ///
      /// @return \f$ p_{Y_{0:N}}(y_{0:N}) \f$
      ///
      double Likelihood_Get(void);

      /// Save the estimation \f$ \{ \hat{X}_{k|k} ,k=0,...N \} \f$
      ///
      /// @param filename 
      ///
      /// @return 0 if everything is ok
      ///
      virtual int Save_X(const char *filename);
protected :

      /// Specific update for each filter
      ///
      /// @param Y The observed sample
      ///
      /// @return 0 if no problem
      ///
      virtual int _update(const dcovector & Y) =0;

      /// Specific init for each filter
      ///
      /// @param Y The observed sample
      ///
      /// @return 0 if no problem
      ///
      virtual int _init(void) =0;

};


///  @brief Abstract class of Gaussian Approximation filters
/// 
///  For  discrete-discrete  models (Gaussian_Nonlinear_Model),  these
///  filters  approximate   the  probability  density   of  the  state
///  transition  \f$ p_{X_k|X_{k-1}}  \f$ and  the probability  of the
///  observation  \f$  p_{Y_k|X_k}  \f$  by gaussian  densities.   The
///  approximation   is   exact   in   the  case   of   linear   model
///  (Gaussian_Linear_Model) and lead  to the discrete-discrete Kalman
///  Filter     (DD_Kalman).     For    other     non-linear    models
///  (Gaussian_Nonlinear_Model)  UKF (Unscented_Kalman_Filter)  or EKF
///  (Extended_Kalman_Filter) can be used.

class GA_Filter :public Filter
{
public :
      /// The current mean \f$ \hat{X}_{k|k}=E[X_k |Y_{0:k}] \f$
      dcovector M;
      /// The current covariance \f$ \hat{P}_{k|k}=E[(X_k-\hat{X}_{k|k})(X_k-\hat{X}_{k|k})] \f$
      dgematrix R;
      /// The prediction \f$ \hat{X}_{k-1|k}=E[X_{k-1} |Y_{0:k}] \f$
      dcovector Xp;
      /// The  prediction covariance \f$ \hat{P}_{k-1|k}=E[(X_k-\hat{X}_{k-1|k})(X_k-\hat{X}_{k-1|k}) ] \f$
      dgematrix Rp;
public :
      /// A constructor
      GA_Filter(void);


      /// A constructor
      ///
      /// @param m A discrete-discrete gaussian non-linear model
      ///
      GA_Filter(Gaussian_Nonlinear_Model *m);
      
      /// Get the current estimation \f$ \hat{X}_{k|k} \f$
      ///
      ///
      /// @return \f$ \hat{X}_{k|k} \f$
      ///
      dcovector Expected_Get(void);
protected :
      virtual int _init(void);
};

/// @brief Abstract class of continuous-discrete filters
/// 
///  For   continusous-discrete   models  (Continuous_Discrete_Model),
///  these filters  approximate the  probability density of  the state
///  transition \f$  p_{X(t_k)|X(t_{k-1})} \f$ and  the probability of
///  the  observation \f$  p_{Y_k|X(t_k)} \f$  by  gaussian densities.
///  The   approximation   is   exact    in   the   case   of   linear
///  continous-discrete  models  (Linear_CD_Model)  and  lead  to  the
///  continuous-discrete   Kalman  Filter   (CD_Kalman).    For  other
///  non-linear models (Continuous_Discrete_Model) Local linearization
///  filter    (LL_Filter)    or    continous-discrete   Filter    EKF
///  (CD_Extended_Kalman_Filter) can be used.

class CD_Filter : public Filter
{
public :
      /// The current mean \f$ \hat{X}_{k|k}=E[X_k |Y_{0:k}] \f$
      dcovector M;
      /// The current covariance \f$ \hat{P}_{k|k}=E[(X_k-\hat{X}_{k|k})(X_k-\hat{X}_{k|k})] \f$
      dgematrix R;
      /// The prediction \f$ \hat{X}_{k-1|k}=E[X_{k-1} |Y_{0:k}] \f$
      dcovector Xp;
      /// The  prediction covariance \f$ \hat{P}_{k-1|k}=E[(X_k-\hat{X}_{k-1|k})(X_k-\hat{X}_{k-1|k}) ] \f$
      dgematrix Rp;
public :
      /// A constructor
      CD_Filter(void);

      /// A constructor
      ///
      /// @param m A discrete-discrete gaussian non-linear model
      ///
      CD_Filter(Continuous_Discrete_Model *m);

      int Save_X(const char *filename);

      /// Get the current estimation \f$ \hat{X}_{k|k} \f$
      ///
      ///
      /// @return \f$ \hat{X}_{k|k} \f$
      ///
      dcovector Expected_Get(void);
protected :
      int _init(void);

};

/// @brief The continuous-discrete kalman filter
/// 
/// Give an exact solution of \f$ \hat{X}_{k|k}\f$ and 
///  \f$ \hat{P}_{k|k} \f$ for continuous-discrete linear models
/// (Linear_CD_Model). 
class CD_Kalman : public CD_Filter
{
      
public :
      CD_Kalman(void);
      /// A constructor
      ///
      /// @param m The continuous discrete model
      ///
      ///
      CD_Kalman(Linear_CD_Model *m);

protected :
      int _update(const dcovector & Y);
};

/// @brief The discrete-discrete kalman filter
/// 
/// Give an exact solution of \f$ \hat{X}_{k|k}\f$ and 
/// \f$ \hat{P}_{k|k} \f$ for discrete-discrete linear models
/// (Gaussian_Linear_Model). 

class DD_Kalman : public GA_Filter
{

public :
      DD_Kalman(void);

      /// A constructor
      ///
      /// @param m The discrete-discrete model
      ///
      DD_Kalman(Gaussian_Linear_Model *m);

protected :
      int _update(const dcovector & Y);
};

#endif
