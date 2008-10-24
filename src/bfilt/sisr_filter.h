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
/// @file   sisr_filter.h
/// @author paul <paul.frogerais@univ-rennes1.fr>
/// @date   Fri Sep 12 18:38:23 2008
/// 
/// @brief   Implementation of Sequential Importance Sampling With Resampling (SISR) Filter.
/// 
/// 
///


#ifndef __SISR__FILTER
#define __SISR__FILTER

#include <bfilt/simulator.h>
#include <bfilt/filter.h>

class Weighted_Sample {

public:
      
      /// @brief The position
      dcovector   Value;
      /// @brief The weight of the sample
      long double Weight;
};


/// @brief the sequential importance sampler used for sisr filter (bootstrap,optimal ...)
/// 
///

class SI_Sampler{
public:
      Simulator *model;
public:

      ///  @brief  The constructor
      ///
      ///
      SI_Sampler(void);
      /**\brief  constructor
       */
      /// @brief The constructor
      ///
      /// @param m A discrete model
      ///
      SI_Sampler(Simulator *m);

      /// @brief draw a set of possible init state
      ///
      /// @param NbSample Number of sample
      ///
      /// @return A set of weighted samples
      ///
      virtual  vector<Weighted_Sample > DrawInitCloud(const int & NbSample) =0;
   
      /// @brief Draw a set of samples from the importance density Xk given Y0:k X0:k-1
      ///
      /// @param Y_k The observation from 0 to k
      /// @param X_km1 The cloud from 0 to km1
      ///
      /// @return A cloud representing the importance density q(Xk|Y0:k,X0:k-1)
      ///
      virtual  vector<Weighted_Sample > Draw(const dcovector & Y_k, 
                                            const vector<Weighted_Sample >  & X_km1) =0;

      /// @brief Modify the weights of cloud for the weighting step in the sisr
      ///
      /// @param cloud The curent coud at k
      /// @param Y_k The observation at k
      /// @param X_km1  the cloud from at km1
      ///
      /// @return The sum of the weights
      ///
      virtual long double Weight(vector<Weighted_Sample > & cloud,
                                 const dcovector  &  Y_k,
                                 const vector<Weighted_Sample >  & X_km1) =0;

};


/// @brief This sampler use the transition  as importance density.
///
///
class Bootstrap_Sampler:public SI_Sampler
{
public :
      Bootstrap_Sampler(void);

      Bootstrap_Sampler(Simulator *m);

      vector<Weighted_Sample > DrawInitCloud(const int & NbSample);
      
      vector<Weighted_Sample > Draw(const dcovector & Y_k, 
                                    const vector<Weighted_Sample > & X_km1);

      long double Weight(vector<Weighted_Sample > & cloud,
                         const dcovector  &  Y_k,
                         const vector<Weighted_Sample >  & X_k);

};

/// @brief This sampler use the optimal importance density.
///
///
class Optimal_Sampler:public SI_Sampler
{
public :

      Optimal_Sampler(void);

      Optimal_Sampler(Opt_Simulator *m);

      vector<Weighted_Sample > DrawInitCloud(const int & NbSample);

      vector<Weighted_Sample > Draw(const dcovector & Y_k, 
                                    const vector<Weighted_Sample >  & X_km1);

      long double Weight(vector<Weighted_Sample > & cloud,
                         const dcovector &Y_k,
                         const vector<Weighted_Sample >  & X_km1);
};






class SISR_Filter : public Filter
{
protected :
      gsl_rng * r;

      int seed;

public:
      /// @brief The  particle clouds at km1
      vector<Weighted_Sample >  cloud_km1; 

      /// @brief The curent particle cloud 
      vector<Weighted_Sample > cloud; 
  
      /// @brief Number of particle
      int NbSample;

      /// @brief the resampling criterion
      float Rc;

      /// @brief the sampler 
      SI_Sampler *Sys;

      /// @brief A constructor
      ///
      ///
      SISR_Filter(void);               

      /// @brief The destructor
      ///
      ///
      ~SISR_Filter(void);

      /// @brief A constructor
      ///
      /// @param Ns  number of sample
      /// @param s a sampler
      ///
      /// @return 
      ///
      SISR_Filter(const int & Ns, SI_Sampler *s);               
    
      /// @brief A constructor
      ///
      /// @param Ns number of sample
      /// @param rc The resampling criterion
      /// @param seed The seed
      /// @param s  a sampler
      ///
      /// @return 
      ///

      SISR_Filter(const int & Ns, const double &rc, const int &seed, SI_Sampler *s);    
      
           
      /// @brief Set the seed of the random number generator of the discret pdf
      ///
      /// @param s The seed 
      ///
      void  SetSeed(const int &s);



      /// @brief The resampling step
      ///
      /// @param Ns 
      ///
      void Resampling(const int &Ns);


      /// \brief get The current cloud
      ///
      ///
      vector<Weighted_Sample > CloudGet(void);


      /// Set the Resampling Criterion
      ///
      /// @param rc The resampling Criterion
      ///
      void SetRc(const float & rc);
  
      dcovector Expected_Get(void);

protected :
      int _update(const dcovector & Yk);

      /// @brief to initialized the first particle cloud of p(X0)
      ///
      ///
      int _init(void);
};

class Bootstrap_Filter : public SISR_Filter
{
      Simulator *sim;
public :
      Bootstrap_Filter(void);
      ~Bootstrap_Filter(void);
      Bootstrap_Filter(const int & Ns,Simulator *s);
      Bootstrap_Filter(const int & Ns,Gaussian_Nonlinear_Model *m);
};

class CD_Bootstrap_Filter : public SISR_Filter
{
      CD_Simulator *sim;
public :
      CD_Bootstrap_Filter(void);
      ~CD_Bootstrap_Filter(void);
      CD_Bootstrap_Filter(const int & Ns,CD_Simulator *s);
      CD_Bootstrap_Filter(const int & Ns,Continuous_Discrete_Model *m);
      CD_Bootstrap_Filter(const int & Ns,Linear_CD_Model *m);
      virtual int Save_X(const char *filename);
};


class OptSISR_Filter : public SISR_Filter
{
      Opt_Simulator *sim;
public :
      OptSISR_Filter(void);
      ~OptSISR_Filter(void);
      OptSISR_Filter(const int & Ns,Opt_Simulator *m);
      OptSISR_Filter(const int & Ns,Gaussian_Nonlinear_Model *m);
};

#endif
