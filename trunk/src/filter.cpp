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

#include <bfilt/filter.h>

Filter::Filter(void)
{
      model=NULL;
      Likelihood=0.;

}

Filter::~Filter(void)
{
}



double Filter::Likelihood_Get(void)
{
      return Likelihood;
}

int Filter::Init(void)
{
      X.clear();

      if(_init())
            {
                  cout<<"BFilt :: Init Filter Problem"<<endl; 
                  return 1;
            }
      else
            {
                  X.push_back(Expected_Get());
                  return 0;
            }
}
int Filter::Update(const dcovector & Y)
{
      if (X.size())
            {
                  if(!_update(Y))
                        {
                              X.push_back(Expected_Get());
                              return 0;
                        }
                  else
                        {     
                              cout<<"BFilt :: Update Filter Problem!"<<endl;
                              return 1;
                        }
                              
            }
      else 
            {
                  cout<<"BFilt :: Filter must be Init"<<endl;
                  return 1;
            }
}

int Filter::Filtering(const vector<dcovector> & Y)
{
      int N = Y.size();
      int k;
      X.resize(N);
      if(model)
            model->Clear();
      if(_init())
            return 1;

      X[0] = Expected_Get();

      if (N)
            {
                  for (k=1; k<N; k++)
                        {
                              if(model)
                                    model->Update();
            
                              if(_update(Y[k]))
                                    return 1;
                              
                              X[k]=Expected_Get();
                        }
                  return 0;
            }
      else 
            {
                  return 1;
            }
}

int Filter::Save_X(const char *filename)
{
      ofstream file(filename);
      int i,j;
      int N=X.size();
      
      if(file)
            {
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
		  cout<<"BFilt :: Pb to access file"<<endl;
		  return 1;
	    }
      return 0;
      
}

GA_Filter::GA_Filter(void)
{
}

GA_Filter::GA_Filter(Gaussian_Nonlinear_Model *m)
{
      model=m;
}

int  GA_Filter::_init(void)
{
      Gaussian_Nonlinear_Model *m=dynamic_cast<Gaussian_Nonlinear_Model *>(model);
      dcovector X0;
      dsymatrix R0;
      m->Get_Init_Parameters(X0,R0);
      M = X0;
      R = R0;
      Likelihood=0.;
      return 0;
}

dcovector GA_Filter::Expected_Get(void)
{
      return M;
}



CD_Filter::CD_Filter(void)
{
}

CD_Filter::CD_Filter(Continuous_Discrete_Model *m)
{
      model=m;
}
int  CD_Filter::_init(void)
{
      Continuous_Discrete_Model *m=dynamic_cast<Continuous_Discrete_Model *>(model);
      dcovector X0;
      dsymatrix R0;
      m->Get_Init_Parameters(X0,R0);
      M = X0;
      R = R0;
      Likelihood=0.;
      return 0;
}

dcovector CD_Filter::Expected_Get(void)
{
      return M;
}

int CD_Filter::Save_X(const char *filename)
{
      Continuous_Discrete_Model *m=dynamic_cast<Continuous_Discrete_Model *>(model);
      ofstream file(filename);
      int i,j;
      int N=X.size();
      double t=0.;
      double Dx=m->Ts;
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
		  cout<<"BFilt :: Pb to access file"<<endl;
		  return 1;
	    }
      return 0;

}



CD_Kalman::CD_Kalman(void)
{
           Likelihood=0.; 
}
CD_Kalman::CD_Kalman(Linear_CD_Model *m):CD_Filter::CD_Filter(m)
{

      dcovector X0;
      dsymatrix R0;
      model=m;
      m->Get_Init_Parameters(X0,R0);

      M = X0;
      R = R0;
}
int CD_Kalman::_update(const dcovector & Y)
{
      Linear_CD_Model *m = dynamic_cast<Linear_CD_Model *>(model);
      dsymatrix Qv=m->Qv;
      const double _2_PI_= 6.283185307;

      dgematrix H=m->H;
      dgematrix Qi;
      dcovector I;
      dgematrix K;

      Xp = m->Get_Mean_Prediction(M);
      Rp = m->Get_Cov_Prediction(R);

      Qi =  H * (Rp * t(H)) + Qv;

      K = Rp* t(H) * i(Qi);
      I= (Y - H*Xp);
      M= Xp + K * I;
      R= Rp - K * H * Rp;

      Likelihood += - 0.5 * I.l * log(_2_PI_) - 0.5 * log(det(Qi)) - 0.5 * CPPL::t(I) * CPPL::i(Qi) * I;

      return 0;

}




DD_Kalman::DD_Kalman(void)
{
           Likelihood=0.; 
}
DD_Kalman::DD_Kalman(Gaussian_Linear_Model *m):GA_Filter::GA_Filter(m)
{
      dcovector X0;
      dsymatrix R0;
      model=m;
      m->Get_Init_Parameters(X0,R0);

      M = X0;
      R = R0;
}
int DD_Kalman::_update(const dcovector & Y)
{
      Gaussian_Linear_Model *m=dynamic_cast<Gaussian_Linear_Model *>(model);

      dsymatrix Qv=m->Qv;
      const double _2_PI_= 6.283185307;

      dgematrix H=m->H;
      dgematrix Qi;
      dcovector I;
      dgematrix K;

      Xp = m->Get_Mean_Prediction(M);
      Rp = m->Get_Cov_Prediction(R);

      Qi =  H * (Rp * t(H)) + Qv;

      K = Rp* t(H) * i(Qi);
      I= (Y - H*Xp);
      M= Xp + K * I;
      R= Rp - K * H * Rp;

      Likelihood += - 0.5 * I.l * log(_2_PI_) - 0.5 * log(det(Qi)) - 0.5 * CPPL::t(I) * CPPL::i(Qi) * I;

      return 0;

}
