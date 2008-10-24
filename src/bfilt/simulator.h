///
/// @file   simulator.h
/// @author paul <paul.frogerais@univ-rennes1.fr>
/// @date   Wed Sep 24 20:34:21 2008
/// 
/// @brief  Implementation of Simulator classes
/// 
/// 
///

#ifndef __SIMULATOR__
#define __SIMULATOR__

#include <bfilt/gaussian_model.h>

class Simulator
{
protected :
      gsl_rng * r;

public :
      Model *model;
      vector<dcovector> X;
      vector<dcovector> Y;
      Simulator(void);
      ~Simulator(void);

      void Set_Seed(const int & s);
      ///  \brief Draw a sample from p(X0)
      ///
      ///
      /// @return A sample from p(X0)
      ///
      virtual  dcovector Draw_Init(void) =0;

      /// @brief Draw a sample from the transition densisty p(Xk|Xk-1)
      ///
      /// @param Xkm1 X(k-1) the preceding state
      ///
      /// @return 
      ///
      virtual  dcovector Draw_Transition(const dcovector & Xkm1) =0;

      /// @brief Calculate the value of the density of probability of Y given X : p(Y|X)
      ///
      /// @param Xk  The state at k
      ///
      /// @return The simulated observation
      ///
      virtual  dcovector Draw_Observation(const dcovector & Xk) =0;

      /// @brief calculate the value of the density of probability of Y given X : p(Y|X)
      ///
      /// @param Y The osbervation
      /// @param X The state
      ///
      /// @return The value of the density
      ///
      virtual  long double Observation_Density(const dcovector & Y, 
                                              const dcovector & X) =0;
  
      /// @brief simulate the markovian model
      ///
      /// @param N The duration
      /// @param X The state trajectory
      /// @param Y The output
      ///
      virtual void  Simulate(const int & N);

      /// @brief Update the simulation of the markovian model
      ///
      /// @param N The duration
      /// @param X The state trajectory
      /// @param Y The output
      ///
      virtual void  Update(void);

      /// @brief Save the simulated state trajectory in filename
      ///
      /// @param filename  The file
      ///
      /// @return 0 if it's ok
      ///
      virtual int Save_X(const char *filename);

      /// @brief Save the simulated observation trajectory in filename
      ///
      /// @param filename The file
      ///
      /// @return  0 if it's ok
      ///
      virtual int Save_Y(const char *filename);

      /// @brief Clear the simulated trajectory X and Y
      ///
      ///      
      void Clear(void);

      /// @brief A pointer for stochastic input
      dcovector (*b)(void *p, gsl_rng * rng); 

protected :
      virtual void _update(void);

      
};
class Opt_Simulator : public Simulator
{
public :
      Opt_Simulator(void);

      /// @brief Draw a sample from the optimal densisty p(Xk|Yk,Xk-1)
      ///
      /// @param Yk The obseration at k
      /// @param Xkm1 X(k-1) the  state value at k-1
      ///
      /// @return A sample from the optimal importance density
      /// 
      virtual  dcovector Draw_Optimal(const dcovector & Yk, const dcovector & Xkm1) =0;

      /// @brief calculate the value of the density of probability of Yk given Xk-1 : p(Yk|Xk-1)
      ///
      /// @param Yk  the osbervation at k
      /// @param Xkm1 The state at k-1
      ///
      /// @return The value of the density p(Yk|Xk-1)
      ///
      virtual  long double Obs_Optimal_Density(const dcovector & Yk, const dcovector &Xkm1) =0;

};
class G_Simulator : public Opt_Simulator
{
public :
      G_Simulator(void);
      G_Simulator(Gaussian_Nonlinear_Model *m);
      dcovector Draw_Init(void);
      dcovector Draw_Transition(const dcovector & Xkm1);
      dcovector Draw_Observation(const dcovector & Xk);
      long double Observation_Density(const dcovector & Y, const dcovector & X);
      dcovector Draw_Optimal(const dcovector & Yk, const dcovector & Xkm1);
      long double Obs_Optimal_Density(const dcovector & Yk, const dcovector &Xkm1);
};

class G_Simulator_WT : public G_Simulator
{
      vector<dcovector> Xt;
      int NB;
      int N;
public :
      G_Simulator_WT(void);
      G_Simulator_WT(Gaussian_Nonlinear_Model *m,const int &NB, const int &N);
      dcovector Draw_Init(void);
};


class CD_Simulator : public Simulator
{
protected :
      int scheme;
      int _a;
public :
      double Dy;
      double Dx;
      
      CD_Simulator(void);
      CD_Simulator(Continuous_Discrete_Model *cd_m,const int &scheme=SRK4, const int &apha=10);

      virtual dcovector Draw_Init(void);
      dcovector Draw_Transition(const dcovector & Xkm1);
      dcovector Draw_Observation(const dcovector & Xk);
      long double Observation_Density(const dcovector & Y, const dcovector & X);
      int Save_X(const char *filename);
      int Save_Y(const char *filename);
      int Simulate(const double & T);

      void Set_Alpha(const int & al);
protected :
      virtual dcovector draw_state(const dcovector &X);
      void  _update(void);
      
};
class CD_Simulator_WT : public CD_Simulator
{
      vector<dcovector> Xt;
      double TB;
      double T;
public :

      CD_Simulator_WT(void);
      CD_Simulator_WT(Continuous_Discrete_Model *cd_m ,const int &scheme,const int &apha, const double &tb,const double &t);
      
      dcovector Draw_Init(void);

};


class LTI_CD_Simulator : public CD_Simulator
{
public :

      LTI_CD_Simulator(void);
      LTI_CD_Simulator(Linear_CD_Model *cd_m, const int &apha=1);

protected :
      dcovector draw_state(const dcovector &X);
};

class LTI_CD_Simulator_WT : public LTI_CD_Simulator
{
      vector<dcovector> Xt;
      double TB;
      double T;
public :

      LTI_CD_Simulator_WT(void);
      LTI_CD_Simulator_WT(Linear_CD_Model *cd_m, const int &apha, const double &tb,const double &t);
      
      dcovector Draw_Init(void);

};



#endif
