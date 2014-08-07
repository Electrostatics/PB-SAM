#ifndef _LAUTIL_H_
#define _LAUTIL_H_

// Functions defined in file "lautil.cpp"
void applyMat(const double * A, const double * X, double * Y, 
	      const double alpha, const double beta, int ma, int na);
void applyMMat(const double * A, const double * X, double * Y, 
	    const double alpha, const double beta, int ma, int nc, int na);
//void solve(const double * A, const double * B, double *X, int ma, int na, int nb);
//double* solve2(const double * A, const double * B, int na, int nb);
//double* computeInvMat(const double * A, int na);

#endif
