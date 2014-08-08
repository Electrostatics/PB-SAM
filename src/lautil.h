#ifndef _LAUTIL_H_
#define _LAUTIL_H_

//! applyMat function
/*!  Function to compute matrix multiplication:
 computes y = alpha*A*x + beta*y
 \param A matrix for multiplication
 \param X vector for multiplication
 \param Y vector for multiplication
 \param alpha a floating point constant for multiplication
 \param beta a floating point constant for multiplication
 \param ma an int for the number of columns in a
 \param na an int for the number of rows in a 	*/
void applyMat(const double * A, const double * X, double * Y, 
	      const double alpha, const double beta, int ma, int na);

//! applyMMat function
/*!  Function to compute matrix multiplication:
 computes y = alpha*A*X + beta*Y
 \param A matrix for multiplication
 \param X matrix for multiplication
 \param Y matrix for multiplication
 \param alpha a floating point constant for multiplication
 \param beta a floating point constant for multiplication
 \param ma an int for the number of columns in a
 \param nc an int for the number of rows in Y
 \param na an int for the number of rows in a 	*/
void applyMMat(const double * A, const double * X, double * Y, 
	    const double alpha, const double beta, int ma, int nc, int na);

#endif
