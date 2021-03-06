#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#include "expansion.h"
#include "triexpan.h"

using namespace::std;

///////////////////////////////////////////////
//!  The CTransform class
/*!	The class that contains all details for all expansion transforms
 The transform is an pm^2 x pn^2 matrix transforming an
 expansion of maximum order pn to an
 expansion of maximum order pm  */
class CTransform
{
public:
 	//!  The CTransform constructor
	/*!	Creates a transform object with given inputs
	 \param mat a pointer to a matrix that performs transforms
	 \param pm an int of the original number of poles
	 \param pn an int of the new desired number of poles
	 \param scale a floating point of the scale of the transform */
	CTransform(const double * mat, int pm, int pn, double scale) :
	m_pm(pm), m_pn(pn), m_scale(scale)
	{
		int mm = pm*pm, nn = pn*pn;
		
		m_mat.resize(mm, vector<double>(nn));
		for (int j = 0; j < nn; j++)
		{
			int k = j*mm;
			for (int i = 0; i < mm; i++, k++)
	    {
	      m_mat[i][j] = mat[k];
	    }
		}
	}
	//!  The CTransform outputRow function
	/*!	Outputs the mnth row of the transform
	 \param mn an int of the desired row number */
  void outputRow(int mm) const
	{
		assert(mm < m_pm*m_pm);
		for (int i = 0; i < m_pn*m_pn; i++)
			cout << m_mat[mm][i] << " ";
		cout << endl;
	}
	//!  The CTransform outputColumn function
	/*!	Outputs the nnth col of the transform
	 \param nn an int of the desired col number */
  void outputColumn(int nn) const
	{
		assert(nn < m_pn*m_pn);
		for (int i = 0; i < m_pm*m_pm; i++)
			cout << m_mat[i][nn] << " ";
		cout << endl;
	}
	//!  The CTransform output function
	/*!	Outputs thetransform */
  void output() const
	{
		for (int i = 0; i < m_pm*m_pm; i++)
			outputRow(i);
	}
	
	
  int get_pm() const
	{return m_pm; }
	
  int get_pn() const
	{return m_pn; }
	
	//!  The CTransform getColExpan function
	/*!	Outputs the an expansion object of the nnth col of
	 the transform
	 \param nn an int of the desired col number */
  CExpan getColExpan(int nn) const
	{
		assert(nn < m_pn*m_pn);
		
		vector<REAL> M;
		M.resize(m_pm*m_pm);
		
		for (int i = 0; i < m_pm*m_pm; i++)
			M[i]= m_mat[i][nn];
		
		CExpan E(M, CRange(m_pm), 1.0);
		return E;
	}
	//!  The CTransform setMatrix function
	/*!	Resets the current transform with given inputs */
  void setMatrix(const double * mat, int pm, int pn)
	{
		assert(pm == m_pm && pn == m_pn );
		int mm = m_pm*m_pm, nn = m_pn*m_pn;
		
		for (int j = 0; j < nn; j++)
		{
			int k = j*mm;
			for (int i = 0; i < mm; i++, k++)
	      m_mat[i][j] = mat[k];		}
		
	}
	
  
	//!  The CTransform apply function
	/*!	Applies a transform to Ein to create Eout  */
  void apply(const CExpan & Ein, CExpan & Eout) const
	{
		const CRange & R = Ein.getRange();
		int p = R.p2();
		assert( p <= m_pn && p <= m_pm);
		Eout.reset(R);
		
		int mm = p*p, nn = p*p;
		
		double s = 1.0;
		int offset = R.p1()*R.p1();
		for (int j = 0; j < nn - offset; j++)
		{
			for (int i = 0; i < mm; i++)
	    {
	      Eout[i] += s*m_mat[i][j+offset]*Ein[j];
				
	    }
		}
	}
	//!  The CTransform apply function
	/*!	Applies a transform to a tricoeff object */
  void apply(const CTriExpan & Ein, CTriExpan & Eout) const
	{
		apply(Ein[0], Eout[0]);
		apply(Ein[1], Eout[1]);
		apply(Ein[2], Eout[2]);
	}
	
  friend ostream& operator<<(ostream& out, const CTransform T);
	
protected:
  double m_scale;
  int m_pm, m_pn;
  vector<vector<double> > m_mat;
};

//!  The CTransform << operator
/*!	Prints out a transform */
inline
ostream& operator<<(ostream& out, const CTransform T)
{
	for (int i = 0; i < T.m_pm*T.m_pm; i++)
	{
	  for (int j = 0; j < T.m_pn*T.m_pn; j++)
		{
			REAL val = T.m_mat[i][j];
			//out<<val<<" ";
			
		}
		//out<<endl;
	}
}

///////////////////////////////////////////////
//!  The CLtoMTransformFull class
/*!	Converts local to multipole transform */
class CLtoMTransformFull : public CTransform
{
	// Full matrix local-to-multipole transform
public:
	CLtoMTransformFull(const double * mat, int pm, int pn, double scale) :
	CTransform(mat, pm, pn, scale) {}	
	void apply(const CLocalExpan & Lin, CMulExpan & Mout) const
	{ CTransform::apply(Lin, Mout); }	
	// Differentiate LTri and MTri later?
	void apply(const CTriExpan & Ein, CTriExpan & Eout) const
	{ CTransform::apply(Ein, Eout); }	
};

///////////////////////////////////////////////
//!  The CMtoLTransformFull class
/*!	Converts multipole to local transform */
class CMtoLTransformFull : public CTransform
{
	// Full matrix multipole-to-local transform
public:
	CMtoLTransformFull(const double * mat, int pm, int pn, double scale) :
	CTransform(mat, pm, pn, scale) {}
	
	void apply(const CMulExpan & Min, CLocalExpan & Lout) const
	{ CTransform::apply(Min, Lout); }
};

///////////////////////////////////////////////
//!  The CLtoLTransformFull class
/*!	Converts Full matrix local-to-local transform */
class CLtoLTransformFull : public CTransform
{
	// Full matrix local-to-local transform
public:
	CLtoLTransformFull(const double * mat, int pm, int pn, double scale) :
	CTransform(mat, pm, pn, scale) {}
	
	void apply(const CLocalExpan & Lin, CLocalExpan & Lout) const
	{ CTransform::apply(Lin, Lout); }
	// Differentiate LTri and MTri later?
	void apply(const CTriExpan & Ein, CTriExpan & Eout) const
	{ CTransform::apply(Ein, Eout); }
};

///////////////////////////////////////////////
//!  The CMtoMTransformFull class
/*!	Converts Full matrix multipole-to-multipole transform */
class CMtoMTransformFull : public CTransform
{
	// Full matrix multipole-to-multipole transform
public:
	CMtoMTransformFull(const double * mat, int pm, int pn, double scale) :
	CTransform(mat, pm, pn, scale) {}
	
	void apply(const CMulExpan & Min, CMulExpan & Mout) const
	{ CTransform::apply(Min, Mout); }	
};


#if 0
inline
CTransform::CTransform(const double * mat, int pm, int pn, double scale) :
m_pm(pm), m_pn(pn), m_scale(scale)
{
	int mm = pm*pm, nn = pn*pn;
	m_mat.resize(mm, vector<double>(nn));
	
	for (int j = 0; j < nn; j++)
	{
		int k = j*mm;
		for (int i = 0; i < mm; i++, k++)
			m_mat[i][j] = mat[k];
	}
}


inline void
CTransform::outputRow(int mm) const
{
	assert(mm < m_pm*m_pm);
	for (int i = 0; i < m_pn*m_pn; i++)
		cout << m_mat[mm][i] << " ";
	cout << endl;
}


inline void
CTransform::outputColumn(int nn) const
{
	assert(nn < m_pn*m_pn);
	for (int i = 0; i < m_pm*m_pm; i++)
		cout << m_mat[i][nn] << " ";
	cout << endl;
}


inline void
CTransform::output() const
{
	for (int i = 0; i < m_pm*m_pm; i++)
		outputRow(i);
}

inline void
CTransform::apply(const CExpan & Ein, CExpan & Eout) const
{
	
	const CRange & R = Ein.getRange();
	int pn = (R.p2() < m_pn ? R.p2() : m_pn);
	Eout.reset(CRange(0, m_pm));
	
	int mm = m_pm*m_pm, nn = pn*pn;
	
	double s = 1.0;
	int offset = R.p1()*R.p1();
	for (int j = 0; j < nn - offset; j++)
	{
		for (int i = 0; i < mm; i++)
		{
			Eout[i] += s*m_mat[i][j+offset]*Ein[j];
			
		}
	}
}

inline void
CTransform::applyTransposed(const CExpan & Ein, CExpan & Eout) const
{
	const CRange & R = Ein.getRange();
	int pm = (R.p2() < m_pm ? R.p2() : m_pm);
	Eout.reset(CRange(0, m_pn));
	
	int mm = pm*pm, nn = m_pn*m_pn;
	
	double s = 1.0;
	int offset = R.p1()*R.p1();
	for (int j = 0; j < mm - offset; j++)
	{
		for (int i = 0; i < nn; i++)
			Eout[i] += s*m_mat[j+offset][i]*Ein[j];
		
		s *= (m_scale * Ein.getScale());
	}
}

#endif


#endif
