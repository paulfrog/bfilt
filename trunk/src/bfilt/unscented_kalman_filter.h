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
/// @file   unscented_kalman_filter.h
/// @author paul <paul.frogerais@univ-rennes1.fr>
/// @date   Fri Sep 12 18:40:15 2008
/// 
/// @brief  Implementation of the unscented kalman filter (See  Julier Uhlmann 1997)
/// 
/// 
///

#ifndef __UKF__FILTER
#define __UKF__FILTER

#include <bfilt/filter.h>


/// \brief The Discrete Unscented Kalman Filter (UKF)
/// 
///

class Unscented_Kalman_Filter : public GA_Filter
{
      /// @brief The square root matrix (cholesky) of Qw
      dgematrix sqrt_Qw;

      /// @brief The square root matrix (cholesky) of Qv
      dgematrix sqrt_Qv;

      /// @brief The sigma points for the state X
      vector<dcovector> sX;

      /// @brief The sigma points for the state noise W
      vector<dcovector> sW;

      /// @brief The sigma points for the observation
      vector<dcovector> sY;

      /// @brief The first weight to compute the mean  
      double w_0;

      /// @brief The first weight to compute the covariance
      double w_0c;

      /// @brief Other weights
      double w;    

public :

      /// @brief A scaled parameter
      float alpha ;


      /// @brief  The constructor
      ///
      /// @param model A gaussian non linear model
      ///
      ///      
      Unscented_Kalman_Filter(Gaussian_Nonlinear_Model * model);

      ///  @brief  The constructor
      ///
      /// @param model A gaussian non linear model
      /// @param a Dispesrion of sigma point
      ///
      ///
      Unscented_Kalman_Filter(Gaussian_Nonlinear_Model * model, const double & a);


protected :
      /// @brief Initialize the sigma points at each update step
      ///
      int SP_Init(void);


      /// @brief Calculate the covaraince between two sets of sigma points
      ///
      /// @param sP1 The first set of sigma point 
      /// @param m1 The mean of the sigma point 
      /// @param sP2 The second set of sigma point 
      /// @param m2 The mean of the second set of sigma point 
      /// @param cov Return the empirical covariance matrix between two sets
      ///
      /// @return 0 if dimensions are ok
      int U_Cov(const vector<dcovector> &sP1, const dcovector &m1,
		const vector<dcovector> &sP2, const dcovector &m2,
		dgematrix & cov);

      /// @brief Calculate the mean of a set of sigma points
      ///
      /// @param sP a set of sigma point
      /// @param mean Return the mean
      ///
      /// @return 0 if dimensions are ok
      ///      
      int U_Mean(const vector<dcovector> &sP, 
		 dcovector & mean);
      

      int _update(const dcovector &Y);

      /// @brief Itialization of the UKF 
      ///
      int _init(void);

};
/// @brief Get a column from a matrix
///
/// @param M The matrix
/// @param k the number of the column
///
/// @return The column vector
///
dcovector get_column(const dgematrix &M, 
		     const int & k);
#endif
