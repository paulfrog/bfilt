///
/// @file   filter_tools.h
/// @author paul <paul.frogerais@univ-rennes1.fr>
/// @date   Fri Sep 12 18:34:54 2008
/// 
/// @brief  Some usefull functions
/// 
/// 
///

#ifndef TOOLS__INCLUDE__FILE
#define TOOLS__INCLUDE__FILE

#include <cpplapack/cpplapack.h>
#include <vector>
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 


using namespace std;
using namespace CPPL;

enum{EULER, HEUN, SRK4, OZAKI, RK4, THGL,RK4_FM};

/// Calculate a eingen decomposition A = P * D * P^-1
///
/// @param A 
/// @param P The matrix of eigen vectors
/// @param D The digonal matrix of eigen values
///
/// @return 0 if ok
///
int eigen_decomposition(const dgematrix &A, dgematrix & P, dgematrix &D); 

/// calulate a square root matrix of A
///
/// @param A a matrix
/// @param sqrtA the square root of A
///
/// @return 0 if ok
///
int sqrtm(const dsymatrix &A, dgematrix & sqrtA);

/// calulate a square root matrix of A
///
/// @param A a matrix
/// @param sqrtA the square root of A
///
/// @return 0 if ok
///
int sqrtm(const dgematrix &A, dgematrix & sqrtA);

/// Compute the exponential of a matrix
///
/// @param A the matrix
///
/// @return the exponential of A
///
dgematrix expm(const dgematrix &A);
double norm_1(const dgematrix &A);
dcovector getPadeCoefficients(const int & m);
dgematrix PadeApproximantOfDegree(const dgematrix A, const int &m);

/** \brief compute the cholesky decompisition
 * \param A the symetric and semi difinite positive matrix
 * \param L the lower matrix
 * \return 1 if A is not semi difinite positive.
 */
int  cholesky(const dsymatrix & A, dgematrix & L);

/** \brief draw from a multivariate gaussian distribution
 * \param X the random result 
 * \param r a pointer to a GSL random number generator
 * \param Q the covariance matrix
 * \return 1 if Q is not semi positive definite
 */
int multivariate_normal_draw(dcovector & X,gsl_rng * r,const dsymatrix & Q);

/** \brief draw from a multivariate gaussian distribution
 * \param X the random result 
 * \param r a pointer to a GSL random number generator
 * \param Q the covariance matrix
 * \return 1 if Q is not semi positive definite
 */
int multivariate_normal_draw(dcovector & X,gsl_rng * r,const dgematrix & Q);

/** \brief return the value of the multinomial gaussian probability density
 * \param X the centered variable
 * \param Q the covariance matrix
 */
long double multivariate_normal_evaluate(const dcovector & X,const  dsymatrix &Q);
/** \brief return the det of a symetric and semi definite posive matrix
 * \param A the matrix
 * \return the det
 */
long double multivariate_normal_evaluate(const dcovector & X,const  dgematrix &Q);
/** \brief return the det of a symetric and semi definite posive matrix
 * \param A the matrix
 * \return the det
 */
double det(const dsymatrix & A);
/** \brief return the det of a symetric and semi definite posive matrix
 * \param A the matrix
 * \return the det
 */
double det(const dgematrix & A);
/** \brief convert a lapack vector into a std vector 
 * \param Z the lapack vector
 * \return the std vector
 */
vector<double> lapack_to_std(const dcovector &Z);
/** \brief convert a std vector into a lapack vector 
 * \param X the std vector
 * \return the lapack vector
 */
dcovector std_to_lapack(const vector<double> &X);
/** \brief convert a lapack matrix into a std vector of std vector 
 * \param X the lapack matrix
 * \return the std vector of std vector
 */
vector<vector<double> > lapack_to_std(const dgematrix & M);
/** \brief convert a std vector of std vector into a lapack matrix
 * \param X the std vecor
 * \return the lapack matrix
 */
dgematrix std_to_lapack(const vector<vector<double> > &M);


int save_dcovector(const  dcovector &X, ostream &s);
int save_dgematrix(const  dgematrix &M, ostream &s);
int save_dsymatrix(const  dsymatrix &M, ostream &s);

int load_dcovector(dcovector &X, istream &s);
int load_dgematrix(dgematrix &M, istream &s);
int load_dsymatrix(dsymatrix &M, istream &s);


vector<double> mean_std_vector(const vector<vector<double > > &X);
dsymatrix cov(const vector<dcovector > &X);
dsymatrix cov(const dgematrix &X);
dcovector mean(const dgematrix &X);
dcovector mean(const vector<dcovector > &X);

int save_signal(const vector<dcovector > &Y, const char *filname, const char *comments="#");
int load_signal(vector<dcovector > &Y, const char *filname);

void shake(gsl_rng *r, vector<dcovector > &X);


template<typename T>
ostream &operator<<(ostream &s, const vector<T> &v){
      int i,N=v.size();
      for(i=0;i<N;i++)
	    s<<v[i]<<" ";
      return s;
}

#endif
