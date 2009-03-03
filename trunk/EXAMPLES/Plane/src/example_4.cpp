#include <bfilt/simulator.h>
#include <bfilt/sisr_filter.h>
#include "plane.h"

// This example illustrate performances of particle filter to highly non-linear 
// filter. The promblem here involves a plane whose the trajectory is a brownian
// motion. This aircraft measure the elevation. The measure of this elevation and 
// an elevation map are then used to estimate the position of the plane.
int main(int argc, char **argv)
{
      int k;
      int i;
      int j;
      vector<Weighted_Sample> cloud;

      ofstream file_c("../data/cloud.dat"); // To save the cloud
      ofstream file_s("../data/state.dat"); // To save the state

      Plane       plane("../data/map_2.dat");       // The plane model
      G_Simulator sim(&plane);                    // To simulate a this model
      
      Bootstrap_Filter filter(100000,&plane);       // A bootstrap filter to estimate the position

      sim.Simulate(150);                     // simulation of 250 samples

      sim.Save_Y("../data/output.dat");     // save the output

      // Here Init() and Update methods are used
      // for filtering because we want to get the cloud 
      // at each step and save it in cloud.dat

      filter.Init();                      // Initialization of the boostrap filter
      
      for (k=0; k<sim.Y.size(); k++)
            {

                  cloud=filter.CloudGet();              // The current cloud is return

                  for(i=0; i<cloud.size(); i++)         // and here it is saved
                        {
                              for(j=0; j<cloud[i].Value.l; j++)
                                    file_c<<cloud[i].Value(j)<<" ";
                              file_c<<endl;
                        }

                  for (j=0; j<sim.X[k].l; j++)       // The state is also saved
                        file_s<<sim.X[k](j)<<" ";
                  file_s<<endl<<endl<<endl;
                  file_c<<endl<<endl;

                  if( filter.Update(sim.Y[k]))              // The filter is then update with a new observation
                        filter.Init();
            }
      file_c.close();
      file_s.close();

      return 0;
}
