#include "lautil.h"

/*
#include "acml.h"

void applyMMat( const double * A, const double * X, double * Y, 
		const double alpha, const double beta, int ma, int nc, int na ) 
{
	const int M = ma;
	const int N = nc;
	const int K = na;
	const int lda = M;
	const int ldb = K;
	const int ldc = M;

	char transA = 'N';
	char transB = 'N';
	dgemm( transA, transB, M, N, K, alpha, const_cast<double*>(A ), 
				lda, const_cast<double*>( X ), ldb, beta, Y, ldc);
	return;
}

void applyMat( const double * A, const double * X, double * Y, 
		const double alpha, const double beta, int ma, int na )
{
	const int N = na;
	const int M = ma;
	const int lda = ma;
	const int incX = 1;
	const int incY = 1;

	char transAc = 'N';
	dgemv( transAc, M, N, alpha, const_cast<double*>(A ), 
				lda, const_cast<double*>( X ), incX, beta, Y, incY);
	return;
}
*/

#include "mkl.h"

// computes y = alpha*A*x + beta*y
void applyMat( const double * A, const double * X, double * Y, 
							const double alpha, const double beta, int ma, int na )
{
  const int N = na;
  const int M = ma;
  const int lda = ma;
  const int incX = 1;
  const int incY = 1;
  const CBLAS_ORDER Order = CblasColMajor;
  const CBLAS_TRANSPOSE TransA = CblasNoTrans;
  cblas_dgemv( Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY );
  return;
}

// computes Y = alpha*A*X + beta*Y
void applyMMat( const double * A, const double * X, double * Y, 
							 const double alpha, const double beta, int ma, int nc, int na ) 
{
  const int M = ma;
  const int N = nc;
  const int K = na;
  const int lda = M;
  const int ldb = K;
  const int ldc = M;
  const CBLAS_ORDER Order = CblasColMajor;
  const CBLAS_TRANSPOSE TransA = CblasNoTrans;
  const CBLAS_TRANSPOSE TransB = CblasNoTrans;
	
  cblas_dgemm(  Order, TransA, TransB, M, N, K, alpha, A, lda, X,  ldb, beta, Y, ldc );
  return;
}


