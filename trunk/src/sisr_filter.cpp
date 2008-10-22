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

#include <bfilt/sisr_filter.h>


SI_Sampler::SI_Sampler(void)
{
}

SI_Sampler::SI_Sampler(Simulator *m)
{
      model=m;
}


Bootstrap_Sampler::Bootstrap_Sampler(void)
{
}

Bootstrap_Sampler::Bootstrap_Sampler(Simulator *m)
{
      model=m;
}

vector<Weighted_Sample > Bootstrap_Sampler::DrawInitCloud(const int & NbSample)
{
      vector<Weighted_Sample > cloud(NbSample);
      int i;

      for(i=0;i<NbSample;i++)
            {
                  cloud[i].Value=model->Draw_Init();
                  cloud[i].Weight=1./(long double)NbSample;
            }

      return cloud;
}
vector<Weighted_Sample > Bootstrap_Sampler::Draw(const dcovector &  Yk, const vector<Weighted_Sample > & X_km1)
{
      int i;
      int N = X_km1.size();
      vector<Weighted_Sample > cloud=X_km1;
      
      for(i=0;i<N;i++)
            {
                  cloud[i].Value=model->Draw_Transition(X_km1[i].Value);
            }

      return cloud;
}

long double Bootstrap_Sampler::Weight(vector<Weighted_Sample > & cloud,
                                      const dcovector  &  Yk,
                                      const vector<Weighted_Sample > & X_km1){
      int         N=cloud.size();
      int         i;
      long double sum=0.;

      for(i=0;i<N;i++)
            {
                  cloud[i].Weight *=  this->model->Observation_Density(Yk, cloud[i].Value);
                  sum+=cloud[i].Weight;
            }

      return sum;
}


Optimal_Sampler::Optimal_Sampler(void)
{
}

Optimal_Sampler::Optimal_Sampler(Opt_Simulator *m):SI_Sampler::SI_Sampler(m)
{
}

vector<Weighted_Sample > Optimal_Sampler::DrawInitCloud(const int & NbSample)
{
      vector<Weighted_Sample > cloud(NbSample);
      int i;
      for(i=0;i<NbSample;i++)
            {
                  cloud[i].Value=model->Draw_Init();
                  cloud[i].Weight=1./(long double)NbSample;	    
            }

      return cloud;
}


vector<Weighted_Sample > Optimal_Sampler::Draw(const dcovector  & Yk,
                                               const vector<Weighted_Sample >  & X_km1)
{
      int i;
      int N=X_km1.size();
      vector<Weighted_Sample > cloud=X_km1;
       Opt_Simulator  *m=dynamic_cast< Opt_Simulator* >(this->model);

      for(i=0;i<N;i++)
            {
                  cloud[i].Value=m->Draw_Optimal(Yk,X_km1[i].Value);
            }

      return cloud;
}


long double Optimal_Sampler::Weight(vector<Weighted_Sample > & cloud,
                                   const dcovector  &  Yk,
                                   const vector<Weighted_Sample >  & X_km1)
      
{
      int i;
      int N=cloud.size();
      long double sum=0.;
       Opt_Simulator  *m=dynamic_cast< Opt_Simulator* >(this->model);

      for(i=0;i<N;i++)
	    {
		  cloud[i].Weight *= m->Obs_Optimal_Density(Yk,X_km1[i].Value);
		  sum+=cloud[i].Weight;
	    }

      return sum;
}


SISR_Filter::SISR_Filter(void)
{
      Rc=1.0;  // always resampling
      Likelihood=0.0;
      NbSample=0;
      seed=time(0);
      Sys=NULL;
      gsl_rng_env_setup(); 
      r = gsl_rng_alloc (gsl_rng_taus);
      gsl_rng_set(r,time(0)); 
}

SISR_Filter::SISR_Filter(const int & Ns, const double &rc, const int &se, SI_Sampler* s)
{
      gsl_rng_env_setup(); 
      r = gsl_rng_alloc (gsl_rng_taus);

      Rc=rc;  // always resampling
      Likelihood=0.0;
      NbSample=Ns;
      Sys=s;
      seed=se;
      model=Sys->model->model;
}
               

SISR_Filter::SISR_Filter(const int & Ns, SI_Sampler *s)
{
      gsl_rng_env_setup(); 
      r = gsl_rng_alloc (gsl_rng_taus);

      NbSample=Ns;
      Sys=s;
      seed=time(0);
      model=Sys->model->model;

      Rc=1.0;  // always resampling
}



SISR_Filter::~SISR_Filter(void)
{
      gsl_rng_free(r);
}


int SISR_Filter::_init(void)
{        
      gsl_rng_set(r,seed);   
      cloud = Sys->DrawInitCloud(NbSample);
      cloud_km1 = cloud;
      Likelihood=0.0;
      return 0;
}



int SISR_Filter::_update(const dcovector & Yk)
{
      long double sum;
      
      // mutation :
      cloud=Sys->Draw(Yk,cloud_km1);

      // ponderation
      Likelihood+=log(Sys->Weight(cloud,Yk,cloud_km1)); 
      
      sum=0.;
      for(int i=0;i<NbSample;i++)
            {
                  sum+=cloud[i].Weight;
            }
      if(sum!=0.)
            {
                  for(int i=0;i<NbSample;i++)
                        {
                              cloud[i].Weight=cloud[i].Weight/sum;
                        }
            }
      else
	    return 1;

      sum=0.;
      for(int i=0;i<NbSample;i++) // calculate the effective size
            {  
                  sum+=cloud[i].Weight * cloud[i].Weight;
            }
      
      long double Neff=1/sum;
      if((float)(Neff/(long double)NbSample)<Rc)   // if degenerescence  then redistribution 
	    { 
		  Resampling(NbSample);
	    }
      cloud_km1 = cloud;
      
      return 0;
}


dcovector SISR_Filter::Expected_Get(void )
{
      dcovector E;

      E= cloud[0].Weight * cloud[0].Value;

      for(int j=1;j<NbSample;j++)
            {
		  E= E + cloud[j].Weight * cloud[j].Value;
            }

      return E;
}


vector<Weighted_Sample > SISR_Filter::CloudGet(void)
{
      return cloud;
}


void SISR_Filter::SetRc(const float & rc)
{
      Rc=rc;
}

void  SISR_Filter::SetSeed(const int &s)
{
      seed = s;
}

void SISR_Filter::Resampling(const int &Ns)
{
      int i;
      vector<int> index(Ns);
      double *weights=new double[NbSample];
      gsl_ran_discrete_t * g;

      for(i=0;i<NbSample;i++)
            {
                  weights[i]=(double)cloud[i].Weight;
            }

      g=gsl_ran_discrete_preproc(NbSample,weights);

      for (i=0;i<Ns;i++)
            {
                  index[i]=gsl_ran_discrete(r,g);
            }

      sort(index.begin(),index.end());

      NbSample = Ns;

      cloud.resize(Ns);
      for(i=0;i<Ns;i++)
            {
                  cloud[i].Value=cloud[index[i]].Value;
                  cloud[i].Weight=(1./(long double)NbSample);
            }

      
      gsl_ran_discrete_free(g);
      
      delete weights;
}

Bootstrap_Filter::Bootstrap_Filter(void)
{
}


Bootstrap_Filter::Bootstrap_Filter(const int & Ns,Simulator *s)
{
      NbSample=Ns;
      sim=s;
      Sys = new  Bootstrap_Sampler(sim);
      model=sim->model;
}

Bootstrap_Filter::Bootstrap_Filter(const int & Ns,Gaussian_Nonlinear_Model *m)
{
      NbSample=Ns;
      sim = new G_Simulator(m);
      Sys = new  Bootstrap_Sampler(sim);
}

Bootstrap_Filter::~Bootstrap_Filter(void)
{
      delete sim;
      delete Sys;
}

CD_Bootstrap_Filter::CD_Bootstrap_Filter(void)
{
}
CD_Bootstrap_Filter::~CD_Bootstrap_Filter(void)
{
      delete sim;
      delete Sys;

}
CD_Bootstrap_Filter::CD_Bootstrap_Filter(const int & Ns,CD_Simulator *s)
{
      NbSample=Ns;
      sim=s;
      model=sim->model;
      Sys = new  Bootstrap_Sampler(s);
}
CD_Bootstrap_Filter::CD_Bootstrap_Filter(const int & Ns,Continuous_Discrete_Model *m)
{
      NbSample=Ns;
      model=m;
      sim = new CD_Simulator(m);
      Sys = new  Bootstrap_Sampler(sim);
}

CD_Bootstrap_Filter::CD_Bootstrap_Filter(const int & Ns,Linear_CD_Model *m)
{
      NbSample=Ns;
      model=m;
      sim = new LTI_CD_Simulator(m);
      Sys = new  Bootstrap_Sampler(sim);
}

int CD_Bootstrap_Filter::Save_X(const char *filename)
{
      ofstream file(filename);
      int i,j;
      int N=X.size();
      double t=0.;
      
      CD_Simulator *s=dynamic_cast< CD_Simulator *> (Sys->model);
      double Dx=s->Dy;
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


OptSISR_Filter::OptSISR_Filter(void)
{
}

OptSISR_Filter::OptSISR_Filter(const int & Ns,Opt_Simulator *s)
{
      NbSample=Ns;
      sim =s;
      Sys = new  Optimal_Sampler(s);
}

OptSISR_Filter::OptSISR_Filter(const int & Ns,Gaussian_Nonlinear_Model *m)
{
      NbSample=Ns;
      sim = new G_Simulator(m);
      Sys = new  Optimal_Sampler(sim);
}

OptSISR_Filter::~OptSISR_Filter(void)
{
      delete sim;
      delete Sys;
}
