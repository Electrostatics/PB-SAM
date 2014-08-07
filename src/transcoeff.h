#ifndef _TRANSCOEFF_H_
#define _TRANSCOEFF_H_

#include "expansion.h"

class CTransCoeff
{
 public:
  CTransCoeff(bool bGrad = false);
  ~CTransCoeff()
    {
      for (int i = m_p; i > 1; i--)
	deallocate();
      m_T[0][0].clear();
      m_T[0].clear();
      m_T.clear();

      if (m_bGrad)
	{
	  m_dT[0][0].clear();
	  m_dT[0].clear();
	  m_dT.clear();
	}
    }
  void reset(REAL rho,REAL kappa, int p, REAL scalei, REAL scaleo);
  void reset(REAL rho,REAL kappa, int p);
  void initScale(REAL scalei, REAL scaleo);

  static void initConstants();
  bool checkT(int l,int n, int m);

  void translate(const CMulExpan & Min, CLocalExpan & Mout, bool tpose) const
    { translate(Min, Mout, m_p, tpose); }
  void translate(const CMulExpan & Min, CLocalExpan & Mout, 
		 int p, bool tpose) const;
  
  void dTranslate(const CMulExpan & Min, CLocalExpan & Mout, bool tpose) const
    { dTranslate(Min, Mout, m_p, tpose); }
  
  void dTranslate(const CMulExpan & Min, CLocalExpan & Mout, 
		  int p, bool tpose) const;
  // void incTranslate(const CMCoeff & Min, CMCoeff & Mout, bool tpose);
  
  REAL computeError(const CMulExpan & tM1, const CMulExpan & tM2, int p);
  int getOrder() const
    { return m_p; }
  
  void incOrder();
  void decOrder();
  
  void outputTrans(int p) const;
  void outputdTrans(int p) const;

  void saveUndo();
  void undo();
  
  void exportMat(double ** T);

 private:
  void allocate();
  void deallocate();
  void reallocate(int p);
  void initParams(REAL d, REAL kappa);
  void computeCoeff(); 
  void computeCoeff_(vector<REAL**> & U);
  void computeCoeff_(vector<vector<vector<REAL> > > & U);
  void computeIncCoeff(); 
  //  void computeIncCoeff_(vector<REAL**> & U);

  static REAL & ALPHA(int n, int m) { return m_alpha[n][m]; }
  static REAL & BETA(int n, int m) { return m_beta[n][m]; }
  static REAL & GAMMA(int n, int m) { return m_gamma[n][N_POLES-1+m]; }
  static REAL & DELTA(int n, int m) { return m_delta[n][N_POLES-1+m]; }
  static bool & EVEN(int n) { return m_even[n]; } 

  static REAL m_alpha[N_POLES*2][N_POLES];
  static REAL m_beta[N_POLES*2][N_POLES];
  static REAL m_gamma[N_POLES*2][2*N_POLES-1];
  static REAL m_delta[N_POLES*2][2*N_POLES-1];
  static bool m_even[4*N_POLES];

  
  REAL & TRANS(int l,int n, int m) { return m_T[l][n][m]; }
  REAL & dTRANS(int l, int n, int m) { return m_dT[l][n][m]; }
  REAL TRANS(int l,int n, int m) const { return m_T[l][n][m]; }
  REAL dTRANS(int l, int n, int m) const { return m_dT[l][n][m]; }
  const REAL FACT(int n) const { return m_fact[n]; } 

  vector<vector<vector<REAL> > > m_T, m_dT;
  //  vector<REAL**> m_T, m_TU, m_dT, m_dTU;
  REAL m_exkid, m_ir, m_d, m_kd;
  REAL m_scalei, m_scaleo;
  vector<REAL> m_K, m_KU;  //  REAL m_K[N_POLES*2], m_KU[N_POLES*2];
  bool m_bGrad;
  REAL m_rbu[4];
  int m_p, m_pU;
  vector<REAL> m_fact;

};

/*
inline void
CTransCoeff::incOrder()
{
  allocate();
  m_p++;

  // increase the bessel function order by two
  CLExpan::INCBESSEL(m_K, m_kd);
  CLExpan::INCBESSEL(m_K, m_kd);

  computeIncCoeff();
}


inline void
CTransCoeff::decOrder()
{
  deallocate();
  m_p--;
}

inline void
CTransCoeff::saveUndo()
{
  
  if (m_bGrad)
    {
      for (int l = 0; l < 2*m_p-1; l++)
	for (int n = 0; n < m_p; n++)
	  for (int m = 0; m <= n; m++)
	    {
	      m_TU[l][n][m] = m_T[l][n][m];
	      m_dTU[l][n][m] = m_dT[l][n][m];
	    }
    }
  else
    {
      for (int l = 0; l < 2*m_p-1; l++)
	for (int n = 0; n < m_p; n++)
	  for (int m = 0; m <= n; m++)
	    m_TU[l][n][m] = m_T[l][n][m];
    }

  m_rbu[0] = m_exkid;
  m_rbu[1] = m_ir;
  m_rbu[2] = m_d;
  m_rbu[3] = m_kd;

  for (int i = 0; i < 2*m_p; i++)
    m_KU[i] = m_K[i];

  m_pU = m_p;
}

inline void
CTransCoeff::undo()
{
  
  if (m_bGrad)
    {
      for (int l = 0; l < 2*m_p-1; l++)
	for (int n = 0; n < m_p; n++)
	  for (int m = 0; m <= n; m++)
	    {
	      m_T[l][n][m] = m_TU[l][n][m];
	      m_dT[l][n][m] = m_dTU[l][n][m];
	    }
    }
  else
  
    {
  
      for (int l = 0; l < 2*m_p-1; l++)
	for (int n = 0; n < m_p; n++)
	  for (int m = 0; m <= n; m++)
	    m_T[l][n][m] = m_TU[l][n][m];
    }

  m_exkid = m_rbu[0];
  m_ir = m_rbu[1];
  m_d = m_rbu[2];
  m_kd = m_rbu[3];

  for (int i = 0; i < 2*m_p; i++)
    m_K[i] = m_KU[i];

  m_p = m_pU;
}
*/






#endif
