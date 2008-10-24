/*!
  \page page1 An AR process

  This is an example on the following auto regressive (AR) process : 
  \f[ X_k = 0.8 X_{k-1} + Wk    \f]
  where, \f$ X_{k} \in \mathcal{R} \f$, \f$ W_k \sim \mathcal{N}(0,0.1) \f$ \par

  This state is then observed by the output \f$ Y_k \in \mathcal{R} \f$ :
  \f[ Y_k = X_k + V_k \f]
  where \f$ V_k \sim \mathcal{N}(0,1) \f$
  
  \section sec1 Define the AR model

  First the ar process must be define as a sister class of Gaussian_Linear_Model.
  \include ar_process.h
  The Gaussian Linear Model are implemented in the following form :
  \f[ X_k = F X_{k-1} + f + G W_k \f]
  \f[ Y_k = H X_{k-1} + h + V_k \f]
  The constructor of AR_Process is then :
  \include ar_process.cpp

  \section sec2 The main program
  In the main program, the model will be first simulted with a specific 
  simulator for gaussian model (G_Simulator). Then the simulated output 
  sequence \f$ y_{0:N}\f$ is given to the input of a discrete-discrete 
  kalman filter (DD_Filter) to estimate the state \f$ \hat{X}_{0:k} \f$.
  First, all this objects are declared :

  \dontinclude example_1.cpp
  \skip main
  \until DD_Kalman

  Then 100 samples are simulated :
  \skipline Simulate
  The kalman filter is apply on the output sequence : 
  \line filter
  You can save the simulated sequences :
  \skip _Y
  \until _X

  and the estimated state :
  \skipline  Save_X

  After compileing and execution, with Gnuplot you can plot :
  \code
  plot 'state.dat' w l, 'estimation.dat' w l, 'output.dat'
  \endcode
  To obtain the following graph :
  \image html "ar_process.jpg" 
  \image latex ar_process.eps
  
  \section sec3 The CMakeList.txt
  \include CMakeLists.txt
*/

/*!  
  \page page2 An Ornstien-Uhlenbeck process

  This  example illustrate  how to  use BFilt  for continuous-discrete
  filtering.  Here the  state  is described  by  the following  linear
  stochastic differential equation :
  \f[
  d \left (
  \begin{array}{c}
  x \\
  \dot{x} \\
  \end{array}
  \right )
  =
  \left (
  \begin{array}{cc}
      0 & 1 \\
      -w_0^2 & -\gamma \\
   \end{array}
   \right )

  \left (
  \begin{array}{c}
  x \\
  \dot{x} \\
 \end{array}
 \right )
 dt
 +
 \left(
 \begin{array}{c}
 0 \\
 b \\
 \end{array}
 \right ) dt
 +
 \left (
 \begin{array}{c}
 0 \\
 g \\
 \end{array}
\right ) dW(t)
   \f]

   Where  W(t) is a  Wiener process,  

   \f$ w_0^2=16, \gamma  = 2, b=8, g=2\f$ and  the initials conditions
   \f$ X_0=(0,0) \f$  and \f$ R_0=diag[0,3]\f$.  
   
   The  state \f$  X(t) =  (x,\dot{x})(t)\f$ is  then observed  by the
   output \f$ Y_k \in \mathcal{R} \f$ :  
   \f[ Y_k = x(t_k) + V_k \f] 

   at discrete time \f$ t_k \f$. The sampling period \f$ T_s = t_{k-1}
   - t_k =  0.2s \f$ and \f$  V_k \sim \mathcal{N}(0,0.001)  \f$.  In fact
   only the position is observed. First this model must be define as a
   sister class  of linear  time invariant continuous  discrete models
   (Linear_CD_Model).
   \include ornstein_uhlenbeck.h  
   The Linear_CD_Model are implemented in the following form : 
   \f[ dX(t) = A X(t)dt + Bdt + C dW(t)  \f] 
   \f[  Y_k =  H X(t_k)  + h  + V_k  \f]
   The  constructor of Ornstein_Uhlenbeck_Model is then :
   \include ornstein_uhlenbeck.cpp

  \section sec2 The main program
  In the main program, the model will be first simulted with a specific 
  simulator for Linear_CD_Model (LTI_CD_Simulator). The simulated output 
  sequence \f$ y_{0:N}\f$ is given to the input of the continuous-discrete 
  kalman filter (CD_Filter) to estimate the state trajectory \f$ \hat{X}_{0:k} \f$.
  First, all this objects are declared :

  \dontinclude example_2.cpp
  \skip main
  \until CD_Kalman

  Then 10 second are simulated :
  \skipline Simulate
  The kalman filter is apply on the output sequence : 
  \line filter
  You can save the simulated sequences :
  \skip _Y
  \until _X

  and the estimated state :
  \skipline  Save_X

  After compileing and execution, with Gnuplot you can plot :
  \code
  plot 'state.dat' w l, 'estimation.dat' w l, 'output.dat'
  \endcode
  To obtain the following graph :
  \image html "ornstein.jpg" 
  \image latex ornstein.eps
  
  \section sec3 The CMakeList.txt
  \include CMakeLists.txt

*/


/*!  \example example_3.cpp

  This example  on the Van der Pol  oscillator shows how to  use BFilt for
  non-linear  continuous-discrete  model. Here the van_der_pol class : \par
  van_der_pol.h
  \include van_der_pol.h
  van_der_pol.cpp
  \include van_der_pol.cpp

  Results can be plotted (here with gnuplot):
  \image html  "van_der_pol.jpg"
  \image latex van_der_pol.eps
  
*/
