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

#include <bfilt/filter_tools.h>
double norm_1(const dgematrix &A)
{
      int i,j;
      double N;
      double S;

      N=0;
      for(i=0; i<A.n;i++)
            N+=fabs(A(i,0));
      
      for(j=1;j<A.m;j++)
            {
                  S=0;
                  for(i=0; i<A.n;i++)
                        S+=fabs(A(i,j));
                  if(S>N)
                        N=S;
            }
      
      return N;            
}
dcovector getPadeCoefficients(const int & m)
{
      dcovector C(m+1);
      
      
      switch(m)
            {
            case 3 :
                  C(0) = 120.;
                  C(1) = 60. ;
                  C(2) = 12. ;
                  C(3) = 1.  ;
                  break;
                  
            case 5 :
                  C(0) = 30240.; 
                  C(1) = 15120.;
                  C(2) = 3360. ;
                  C(3) = 420.  ;
                  C(4) = 30.   ;
                  C(5) = 1.    ;
                  break;
            case 7 :
                  C(0) = 17297280.;
                  C(1) = 8648640. ;
                  C(2) = 1995840. ;
                  C(3) = 277200.  ;
                  C(4) = 25200.   ;
                  C(5) = 1512.    ;
                  C(6) = 56.      ;
                  C(7) = 1.       ;
                  
                  break;
            case 9 :
                  C(0) = 17643225600.;
                  C(1) = 8821612800. ;
                  C(2) = 2075673600. ;
                  C(3) = 302702400.  ;
                  C(4) = 30270240.   ;
                  C(5) = 2162160.    ;
                  C(6) = 110880.     ;
                  C(7) = 3960.      ;
                  C(8) = 90.         ;
                  C(9) = 1.          ;
                  break;
            case 13 :
                  C(0) = 64764752532480000. ;
                  C(1) = 32382376266240000. ;
                  C(2) = 7771770303897600.  ;
                  C(3) = 1187353796428800.  ;
                  C(4) = 129060195264000.   ;
                  C(5) = 10559470521600.    ;
                  C(6) = 670442572800.      ;
                  C(7) = 33522128640.       ;
                  C(8) = 1323241920.        ;
                  C(9) = 40840800.          ;
                  C(10) = 960960.           ;
                  C(11) = 16380.            ;
                  C(12) = 182.              ;
                  C(13) = 1.                ;
            default :
                  break;

            }
      return C;
}
dgematrix PadeApproximantOfDegree(const dgematrix A, const int &m)
{
      int  n = A.n;
      dcovector c = getPadeCoefficients(m);
      int N = (m+1)/2;
      int i,j;
      vector<dgematrix> Apower(N);
      dgematrix U(n,n);
      dgematrix V(n,n);
      dgematrix F;
      dgematrix A2,A4,A6,I(n);
      
      if(m<12)
            {
                  Apower[0].resize(n,n);
                  Apower[0].identity();
                  Apower[1] = A * A ;
                  
                  for (i=2; i<N; i++)
                        Apower[i]=Apower[i-1] * Apower[1];
                  
                  U. zero();
                  for(j=m+1; j >= 2; j-=2)
                        U += c(j-1) * Apower[j/2-1];
                  U=A*U;
                  V.zero();
                  for(j=m  ; j >= 1; j-=2)
                        V += c(j-1) * Apower[(j+1)/2-1];
                  
                  F = (U+V)*CPPL::i(-U + V) ;
                  
            }
      else 
            {
                  I.identity();
                  A2 = A  * A ;
                  A4 = A2 * A2;
                  A6 = A2 * A4;
                  U = A * (A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*I );
                  V = A6*(c(12)*A6 + c(10)*A4 + c(8)*A2)+ c(6)*A6 + c(4)*A4 + c(2)*A2 + c(0)*I;
                  
                  F = (U+V)*CPPL::i(-U + V) ;
            }
            
      return F;
}
dgematrix expm(const dgematrix &A)
{
      int  m_vals[5] = {3,5,7,9,13};
      long double theta[5] = {1.495585217958292e-002,
                              2.539398330063230e-001,
                              9.504178996162932e-001,
                              2.097847961257068e+000,
                              5.371920351148152e+000};
      
      
      double normA = norm_1(A);
      int i;
      int s;
      int p;
      dgematrix F;
      if (normA <= theta[4])
            {
                  for (i=0;i<5;i++)
                        if(normA <= theta[i])
                              {
                                    F = PadeApproximantOfDegree(A,m_vals[i]);
                                    break;
                              }
            }
      else
            {
                  s = (int)(log2(normA/theta[4])) + 1;
                  p = (int)pow(2.,s);
                  F = PadeApproximantOfDegree(A/p,m_vals[4]);
                  for( i=0 ; i < s ; i++)
                        F*=F;
                  
                        
            }
      return F;
 }


double det(const dsymatrix & A){
      dgematrix L;
      cholesky(A,L);
      double d=1.;
      int i,N=A.n;
      for(i=0;i<N;i++)
	    d*=L(i,i)*L(i,i);
      return d;
}

double det(const dgematrix & A){
      dgematrix L;
      cholesky(A,L);
      double d=1.;
      int i,N=A.n;
      for(i=0;i<N;i++)
	    d*=L(i,i)*L(i,i);
      return d;
}


int cholesky(const dgematrix & A, dgematrix &L){
      int n,N=A.m;
      int i,j;
      if (A.n != N) {cout << "matrix not square !"<<endl; return 1;}
      dsymatrix S;
      S.resize(N);
      for (j=0;j<N;j++)
	    for(i=j;i<N;i++)
		  S(i,j)=A(i,j);
      return  cholesky(S,L);
}
int cholesky(const dsymatrix & A, dgematrix &L){

      const int N=A.n;
      L.resize(N,N);
      L.zero();
  
      int i,j,k;
      double sum,x;

      x=(A(0,0));
      if(x<=0.){return 1;}
      L(0,0)=sqrt(x);
      for (i=1;i<N;i++)
	    L(i,0)=A(0,i)/L(0,0);
      i=1;
      for (i=1;i<N;i++){
	    sum=0.;
	    for (k=0;k<i;k++)
		  sum+=L(i,k)*L(i,k);
	    x=(A(i,i)-sum);
	    if(x<=0.){ return 1;}
	    L(i,i)=sqrt(x);
	    for (j=i+1;j<N;j++){
		  sum=0.;
		  for (k=0;k<i;k++)
			sum+=L(i,k)*L(j,k);
		  L(j,i)=(A(i,j) - sum)/L(i,i);
	    }
      }

      return 0; 
}

int multivariate_normal_draw(dcovector & X, gsl_rng * r, const dsymatrix & Q){
      X.resize(Q.n);
      dgematrix L;
      int i;
      for (i=0;i<X.l;i++)
	    X(i)=gsl_ran_gaussian(r,1.);
      if(  cholesky(Q,L)){
	    X.zero();
	    return 1;
      }
      X=L*X;
      return 0;
}

int multivariate_normal_draw(dcovector & X, gsl_rng * r, const dgematrix & Q){
      X.resize(Q.n);
      dgematrix L;
      int i;
      for (i=0;i<X.l;i++)
	    X(i)=gsl_ran_gaussian(r,1.);
      if(  cholesky(Q,L)){
	    X.zero();
	    return 1;
      }
      X=L*X;
      return 0;
}

long double multivariate_normal_evaluate(const dcovector & X, const dsymatrix & Q){
      int N=X.l;
      long double x;
      dgematrix  Qinv=i(Q);

      const double _2_PI_= 6.283185307;
      x =t(X)* Qinv  * X;
      long double K = 1./(pow(_2_PI_,N-1)*sqrt(_2_PI_)* sqrt(det(Q)) );
      return K*exp(-0.5*x);
}

long double multivariate_normal_evaluate(const dcovector & X, const dgematrix & Q){
      int N=X.l;
      double x;
      dgematrix  Qinv=i(Q);
      const double _2_PI_= 6.283185307;
      x =t(X)* Qinv  * X;
      double K = 1/(pow(_2_PI_,N-1)*sqrt(_2_PI_)* sqrt(det(Q)) );
      return K*exp(-0.5*x);
}

vector<double> lapack_to_std(const dcovector & Z){
      vector<double> X(Z.l);
      int i; 
      for (i=0;i<Z.l;i++)
	    X[i]=Z(i);
      return X;
}
dcovector std_to_lapack(const vector<double> &X){
      dcovector Z(X.size());
      int i;
      for (i=0;i<Z.l;i++)
	    Z(i)=X[i];
      return Z;
}


vector<vector<double> > lapack_to_std(const dgematrix & M){
      int m=M.m;
      int n=M.n;
      int i,j;
      vector<vector<double> > A(m);
      for (i=0;i<m;i++){
	    A[i].resize(n);
	    for(j=0;j<n;j++)
		  A[i][j]=M(i,j);
      }
      return A;
}

dgematrix std_to_lapack(const vector<vector<double> > &M){
      int n=M.size();
      int m=0;
      int i,j;
      if(n!=0)
	    m=M[0].size();
      dgematrix A(m,n);

      for (i=0;i<m;i++)
	    for(j=0;j<n;j++)
		  A(i,j)=M[j][i];
      return A;
}


int save_dcovector(const  dcovector &X, ostream &s){

      s<<X.l<<endl;

      s<<X;

      return 0;
}


int save_dgematrix(const  dgematrix &M, ostream &s){

      s << M.m <<" "<<M.n<<endl;

      s << M;

      return 0;
}


int save_dsymatrix(const  dsymatrix &M, ostream &s){

      int i,j;

      s <<M.n<<endl;

      for(i=0; i<M.n; i++){

	    for(j=0; j<M.n; j++)

		  s<<M(i,j)<<" " ;
    
	    s<<endl;

      }

      return 0;
}


int load_dcovector(dcovector &X, istream &s){

      int i , N;

      s>>N;

      X = dcovector(N);
  
      for(i=0;i<N;i++)
    
	    s>>X(i);
  
      return 0;
}


int load_dgematrix(dgematrix &A, istream &s){

      int i,j,M,N;;
  
      s>>M;
  
      s>>N;
  
      A=dgematrix(M,N);

      for (i=0; i<M; i++)
    
	    for(j=0; j<N; j++)
    
		  s>>A(i,j);

      return 0;
}
int load_dsymatrix(dsymatrix &A, istream &s){

      int i,j,N;;

      s>>N;

      A = dsymatrix(N);
  
      for (i=0; i<N; i++)
    
	    for(j=0; j<N; j++)
      
		  s>>A(i,j);

      return 0;
}

#ifdef __MPI_USED

MPI_MESSAGE& operator <<(MPI_MESSAGE& mpi_message, const dcovector &X)
{
      int i;

      mpi_message<<X.l;

      for (i=0; i<X.l; i++)
	    mpi_message<<X(i);

      return mpi_message;

}

void operator >>(MPI_MESSAGE& mpi_message, dcovector &X)
{
      int i;
      int N;
      mpi_message>>N;
      X=dcovector(N);
      
      for (i=0; i<N; i++)
	    mpi_message>>X(i);
}


MPI_MESSAGE& operator <<(MPI_MESSAGE& mpi_message, const dgematrix &M)
{
      int i,j;

      mpi_message<<M.m;
      mpi_message<<M.n;
      for (i=0; i<M.m; i++)
	    for(j=0; j<M.n; j++)
		  mpi_message<<M(i,j);
      
      return mpi_message;
}

void operator >>(MPI_MESSAGE& mpi_message, dgematrix &A)
{
      int i,j;
      int M,N;
      mpi_message>>M;
      mpi_message>>N;
      A=dgematrix(M,N);
      for (i=0; i<A.m; i++)
	    for(j=0; j<A.n; j++)
		  mpi_message>>A(i,j);
}

MPI_MESSAGE& operator <<(MPI_MESSAGE& mpi_message, const dsymatrix& M)
{
      int i,j;

      mpi_message<<M.n;
      for (i=0; i<M.n; i++)
	    for(j=0; j<M.n; j++)
		  mpi_message<<M(i,j);
      
      return mpi_message;

}

void operator >>(MPI_MESSAGE& mpi_message, dsymatrix &M)
{

      int i,j;
      int N;
      mpi_message>>N;
      M=dsymatrix(N);
      for (i=0; i<M.n; i++)
	    for(j=0; j<M.n; j++)
		  mpi_message>>M(i,j);
}


#endif


int save_signal(const vector<dcovector > &Y, const char *filename, const char *s){

      ofstream file(filename);
      int i,j;
      int D;
      int N=Y.size();

      file<<s<<endl;

      if(file){
	    for (i=0;i<N;i++)
                  {
                        for (j=0; j<Y[0].l; j++)
                              file<<Y[i](j)<<" ";
                        file<<endl;
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


int load_signal(vector<dcovector > &Y, const char *filename){
      ifstream file(filename);
      Y.clear();
      string line;
      dcovector y(1);
      int g;

      if(file){
	    do{
		  g=file.tellg();
		  getline(file,line);
	    }while(line.find("#")==0);
	    file.seekg(g);
	    while(!file.eof()){
		  file>>y(0);
		  Y.push_back(y);
	    }
      }

      else
	    {
		  cout<<"ERROR : Pb to access file"<<endl;
		  return 1;
	    }
      return 0;
}


vector<double> mean_std_vector(const vector<vector<double > > &X){
      vector<double> MOY;
      int N = X.size();
      int M;
      int i,j;

      if(N){
	    M=X[0].size();
	    MOY.resize(M);
	    for(i=0; i<M; i++)
		  {
			for(j=0; j<N; j++)
			      MOY[i]+=X[j][i];
			MOY[i]/=N;
		  }
      }
      return MOY;
}


void shake(gsl_rng *r, vector<dcovector > &X){
      vector<dcovector > Y=X;
      int i , N= X.size();
      int *k=new int[N];
      for (i=0; i<N; i++)
	    k[i]=i;      

      gsl_ran_shuffle(r,k,N,sizeof(int));
      
      for(i=0;i<N;i++)
	    X[i]=Y[k[i]];

	    
      delete k;
}


dcovector mean(const dgematrix &X){
      dcovector m(X.m);
      int i, j;
      m.zero();
      for (i=0; i<X.m; i++)
	    {
		  for (j=0;j<X.n; j++)
			m(i)+=X(i,j);
		  m(i)/=X.n;
	    }
      return m;
}

dcovector mean(const vector<dcovector > &X)
{
     
      dcovector m;
      int i, j;
      int N=X.size();
      if(N!=0)
            {
                  m.resize(X[0].l);
                  m.zero();
                  for (i=0; i<X[0].l; i++)
                        {
                              for (j=0;j<N; j++)
                                    m(i)+=X[j](i);
                              m(i)/=N;
                        }
            }
      return m;
}

dsymatrix cov(const vector<dcovector > &X)
{
      dgematrix U;
      dsymatrix Q; 
      dcovector m=mean(X);
      int i, j, k;
      int N=X.size();

      if (N!=0)
            {
                  U.resize(X[0].l,N);
                  Q.resize(X[0].l);
                  for (i=0; i<X[0].l; i++)
                        for (j=0;j<N; j++)
                              U(i,j)=X[j](i)-m(i);
                  
                  for (i=0; i<Q.n; i++)
                        for(j=i; j<Q.n; j++)
                              {
                                    Q(i,j)=0;
                                    for (k=0;k<U.n; k++)
                                          Q(i,j) += U(i,k)*U(j,k);
                                    Q(i,j)/=U.n;
                              }
            }
      return Q;
}


dsymatrix cov(const dgematrix &X)
{
      dgematrix U=X;
      dsymatrix Q(X.m); 
      dcovector m=mean(X);
      int i, j, k;
      for (i=0; i<X.m; i++)
	    for (j=0;j<X.n; j++)
		  U(i,j)=X(i,j)-m(i);
      
      for (i=0; i<Q.n; i++)
	    for(j=i; j<Q.n; j++)
		  {
			Q(i,j)=0;
			for (k=0;k<X.n; k++)
			      Q(i,j) += U(i,k)*U(j,k);
			Q(i,j)/=X.n;
		  }
      return Q;
}
