///
/// @file   extended_kalman_filter.h
/// @author paul <paul.frogerais@univ-rennes1.fr>
/// @date   Fri Sep 12 18:34:06 2008
/// 
/// @brief  Implementation of the Extended Kalman Filter for discrete models
/// 
/// 
///
#ifndef _E_KALMAN__H
#define _E_KALMAN__H

#include <bfilt/filter.h>



class Extended_Kalman_Filter : public GA_Filter
{
public:
      Extended_Kalman_Filter(void);
      Extended_Kalman_Filter(Gaussian_Nonlinear_Model *m);

protected :
      int  _update(const dcovector &Y);
};


class CD_Extended_Kalman_Filter : public CD_Filter
{
public:
      int Scheme;
      CD_Extended_Kalman_Filter(void);
      CD_Extended_Kalman_Filter(Continuous_Discrete_Model *m, const int & sh=RK4);

protected :
      int  _update(const dcovector &Y);
      void _thgl__prediction(dcovector &M, dgematrix &P);
      void _euler_prediction(dcovector &M, dgematrix &P);
      void _rk4___prediction(dcovector &M, dgematrix &P);
      void _heun__prediction(dcovector &M, dgematrix &P);
};



#endif
