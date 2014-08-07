#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#include "expansion.h"
#include "triexpan.h"

using namespace::std;

////////////////////////////////////////////////
//    CTransform
///////////////////////////////////////////////
class CTransform
{
  // A base class for all expansion transforms
  // The transform is and pm^2 x pn^2 matrix transforming a an
  // expansion of maximum order pn to an expansion of maximum order
  // pm
 public:
  CTransform(const double * mat, int pm, int pn, double scale) :
    m_pm(pm), m_pn(pn), m_scale(scale)
    {

      int mm = pm*pm, nn = pn*pn;

      m_mat.resize(mm, vector<double>(nn));
      //  cout<<"in transform ... "<<endl;
      for (int j = 0; j < nn; j++)
	{
	  int k = j*mm;
	  for (int i = 0; i < mm; i++, k++)
	    {
	      //      cout<<j<<" "<<i<<" "<<k<<" "<<mat[k]<<endl;
	      m_mat[i][j] = mat[k];
	    }
	}
      //      cout<<"exited CTransform normally."<<endl;
    }

  void outputRow(int mm) const
    {
      assert(mm < m_pm*m_pm);
      for (int i = 0; i < m_pn*m_pn; i++)
	cout << m_mat[mm][i] << " ";
      cout << endl;
    }
  void outputColumn(int nn) const
    {
      assert(nn < m_pn*m_pn);
      for (int i = 0; i < m_pm*m_pm; i++)
	cout << m_mat[i][nn] << " ";
      cout << endl;
    }
  void output() const
    {
      for (int i = 0; i < m_pm*m_pm; i++)
	outputRow(i);
    }


  int get_pm() const
    {return m_pm; }

  int get_pn() const
    {return m_pn; }

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

  void setMatrix(const double * mat, int pm, int pn)
    {
      assert(pm == m_pm && pn == m_pn );
      int mm = m_pm*m_pm, nn = m_pn*m_pn;

      for (int j = 0; j < nn; j++)
	{
	  int k = j*mm;
	  for (int i = 0; i < mm; i++, k++)
	      m_mat[i][j] = mat[k];
	}
      
    }

  
  // apply: output must have same range as input
  void apply(const CExpan & Ein, CExpan & Eout) const
    {
      const CRange & R = Ein.getRange();
      int p = R.p2();
      /*
      if( !(p <= m_pn && p <= m_pm)) 
	cout <<"p "<<p<<" m_pn "<<m_pn<< " m_pm "<< m_pm<<endl;
      */
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
	  //	  s *= (m_scale * Ein.getScale());
	}
    }

  // Itay's version
  /*
  void apply(const CExpan & Ein, CExpan & Eout) const
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
	       //if(fabs(m_mat[i][j+offset])>1e-15)
	//		cout<<i<<" "<<j+offset<<" "<<Ein[j]<<
		//	" "<<m_mat[i][j+offset]<<" "<<Eout[i]<<endl;
	    }
	  //	  s *= (m_scale * Ein.getScale());
	}
    }
*/

  // NOT USED
  /*
  void applyTransposed(const CExpan & Ein, CExpan & Eout) const
    {
      const CRange & R = Ein.getRange();
      int pm = (R.p2() < m_pm ? R.p2() : m_pm);
      Eout.reset(CRange(0, m_pn));
      
      int mm = pm*pm, nn = m_pn*m_pn;
      
      double s = 1.0;
      int offset = R.p1()*R.p1();
      for (int j = 0; j < mm -offset; j++)
	{
	  for (int i = 0; i < nn; i++)
	    Eout[i] += s*m_mat[j+offset][i]*Ein[j];
	  
	  s *= (m_scale * Ein.getScale());
	}
    }
  */
  // By Enghui
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

inline
ostream& operator<<(ostream& out, const CTransform T)
{
      for (int i = 0; i < T.m_pm*T.m_pm; i++)
	{
	  for (int j = 0; j < T.m_pn*T.m_pn; j++)
	    {
	     REAL val = T.m_mat[i][j];
	     out<<val<<" ";

	    }
	  out<<endl;
	}
}

////////////////////////////////////////////////
//    CLtoMTransform
///////////////////////////////////////////////
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

////////////////////////////////////////////////
//    CMtoLTransform
///////////////////////////////////////////////
class CMtoLTransformFull : public CTransform
{
  // Full matrix multipole-to-local transform
 public:
  CMtoLTransformFull(const double * mat, int pm, int pn, double scale) :
    CTransform(mat, pm, pn, scale) {}

  void apply(const CMulExpan & Min, CLocalExpan & Lout) const
    { CTransform::apply(Min, Lout); }
};

////////////////////////////////////////////////
//    CLtoLTransform
///////////////////////////////////////////////
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

////////////////////////////////////////////////
//    CMtoMTransform
///////////////////////////////////////////////
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
      //      s *= (m_scale * Ein.getScale());
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
