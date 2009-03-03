#include "plane.h"

Plane::Plane(const char * filename)
{
      int x,y,z;
      ifstream file(filename);

      if(file)
            {
                              file>>xmax;
                              file>>ymin;
                              file>>z;
                              Map.push_back(z);

                  while(!file.eof())
                        {
                              file>>x;
                              file>>y;
                              file>>z;
                              Map.push_back(z);

                              if(x>xmax)
                                    xmax=x;
                              if(x<xmin)
                                    xmin=x;
                              if(y>ymax)
                                    ymax=y;
                              if(y<ymin)
                                    ymin=y;
                              
                        }
                  file.close();
            }
      else 
            {
                  cout<<"Plane :: error file"<<endl;
            }

      Qw.resize(2);
      Qw.identity();
      Qv.resize(1);
      Qv.identity();
      Qv*=5.;
      R0.resize(4);
      R0.zero();
      R0(0,0)=10.;
      R0(1,1)=10.;
      R0(2,2) = .001;
      R0(3,3) = 0.0001;
      X0.resize(4);
      X0(0)=120.;
      X0(1)=20.;
      X0(2)=1.5;
      X0(3)=2.35;


      sigv = .001;
      sigc = 0.03;
            
}

dcovector Plane::State_Function(const dcovector &X, const dcovector &W)
{
      dcovector U(4);

      U(0) = X(0) + X(2) * cos(X(3));
      U(1) = X(1) + X(2) * sin(X(3));
      U(2) = X(2) + sigv * W(0);
      U(3) = X(3) + sigc * W(1);
      if (U(0)>xmax)
            {
                  U(0)=xmax;
                  U(3)=3.14-U(3);
            }
      if (U(1)>ymax)
            {
                  U(1)=ymax;
                  U(3)=-U(3);
            }
      if (U(0)<xmin)
            {
                  U(0)=xmin;
                  U(3)=3.14-U(3);
            }
      if (U(1)<ymin)
            {
                  U(1)=ymin;
                  U(3)=-U(3);
            }

      return U;
}
dcovector Plane::Observation_Function(const dcovector & X)
{
      int x = (int)(X(0));
      int y = (int)(X(1));
      int j= (xmax+1)*y + x;
      dcovector Y(1);
      Y(0)=Map[j];
      return Y;
}

