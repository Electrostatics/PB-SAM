#include "transcoeff.h"

REAL CTransCoeff::m_alpha[N_POLES*2][N_POLES];
REAL CTransCoeff::m_beta[N_POLES*2][N_POLES];
REAL CTransCoeff::m_gamma[N_POLES*2][2*N_POLES-1];
REAL CTransCoeff::m_delta[N_POLES*2][2*N_POLES-1];
bool CTransCoeff::m_even[4*N_POLES];

/******************************************************************/
/******************************************************************/
/**
 * Initialize the rotation coefficients (includes coeff equs from
 * 2006 paper
 ******************************************************************/
void
CTransCoeff::initConstants()
{
  REAL kapsqr = (CRExpan::KAPPA == 0 ? 0 : 1.0); // i.e. set lambda = 1/kappa
	
  for (int n = 0; n < 2*N_POLES; n++)
    for (int m = 0; m < N_POLES; m++)
		{
			ALPHA(n,m) = sqrt((REAL)(n+m+1)*(n-m+1));						// Given as alpha(n,m) in Lotan, 2006 (eq 1.9)
			BETA(n,m)  = kapsqr*ALPHA(n,m)/((2*n+1)*(2*n+3));		// Given as beta(n,m) in Lotan, 2006 (eq 1.9)
			
			GAMMA(n,m) = sqrt((REAL)(n-m-1)*(n-m));							// Given as eta(n,m) in Lotan, 2006 (eq 1.9)
			DELTA(n,m) = kapsqr*GAMMA(n,m)/((2*n+1)*(2*n-1));		// Given as mu(n,m) in Lotan, 2006 (eq 1.9)
			if (m != 0)
			{
				GAMMA(n,-m) = -sqrt((REAL)(n+m-1)*(n+m));					// Given as eta(n,m) in Lotan, 2006 (eq 1.9)
				DELTA(n,-m) = kapsqr*GAMMA(n,-m)/((2*n+1)*(2*n-1));// Given as mu(n,m) in Lotan, 2006 (eq 1.9)
			}
		}
	
  EVEN(0) = true;
  for (int n = 1; n < 4*N_POLES; n++)
    EVEN(n) = !EVEN(n-1);	
}

/******************************************************************/
/******************************************************************/
/**
 *  		Constructing a transcoeff class, using input of whether or
 not to compute the gradient
 ******************************************************************/
CTransCoeff::CTransCoeff(bool bGrad) : m_bGrad(bGrad), m_p(1)
{
  m_T.resize(1);
  m_T[0].resize(1);
  m_T[0][0].resize(1);
	
  if (m_bGrad)
	{
		m_dT.resize(1);
		m_dT[0].resize(1);
		m_dT[0][0].resize(1);
	}
  
}

/******************************************************************/
/******************************************************************/
/**
 *  Allocating space for the T and dT matrices
 ******************************************************************/
void
CTransCoeff::allocate()
{
  assert(m_T.size() <= 2*N_POLES-3);
	
  int p = m_T.size();
  m_T.resize(p+2);
  m_T[p].resize(p+1);
  m_T[p+1].resize(p+2);
  for (int j = 0; j < p+1; j++)
    m_T[p][j].resize(p+1);
  for (int j = 0; j < p+2; j++)
    m_T[p+1][j].resize(p+2);
	
  if (m_bGrad)
	{
		m_dT.resize(p+2);
		m_dT[p].resize(p+1);
		m_dT[p+1].resize(p+2);
		for (int j = 0; j < p+1; j++)
			m_dT[p][j].resize(p+1);
		for (int j = 0; j < p+2; j++)
			m_dT[p+1][j].resize(p+2);
	}
  
}

/******************************************************************/
/******************************************************************/
/**
 *   Deallocating space for the T and dT matrices
 *****************************************************************/
void
CTransCoeff::deallocate()
{
  assert(m_T.size() > 1);
  int p = m_T.size();
	
  for (int j = 0; j < p; j++)
    m_T[p-1][j].clear();
  for (int j = 0; j < p-1; j++)
    m_T[p-2][j].clear();
  m_T[p-1].clear();
  m_T[p-2].clear();
  m_T.resize(p-2);
	
  if (m_bGrad)
	{
		for (int j = 0; j < p; j++)
			m_dT[p-1][j].clear();
		for (int j = 0; j < p-1; j++)
			m_dT[p-2][j].clear();
		m_dT[p-1].clear();
		m_dT[p-2].clear();
		m_dT.resize(p-2);
	}
}

/******************************************************************/
/******************************************************************/
/**
 *  Function to reallocate space in T and dT matrices
 calls either de or allocate according to m_p
 ******************************************************************/
void
CTransCoeff::reallocate(int p)
{
	
	assert(p <= N_POLES && p >= 1);
  if (p - m_p > 0)
    for (int i = m_p; i < p; i++)
      allocate();
  else
    for (int i = m_p; i > p; i--)
      deallocate();
	
	
}

/******************************************************************/
/******************************************************************/
/**
 *  Initialize parameters in TransCoeff object
 ******************************************************************/
void
CTransCoeff::initParams(REAL d, REAL kappa)
{
  m_d = d;
  m_kd = kappa * m_d;
  
  if (d > 0.0)
	{
		m_ir = 1.0/m_d;
		m_exkid = exp(-m_kd)*m_ir;
	}
	
  m_K.assign(2*m_p, 1.0);
  CLExpan::BESSEL(m_K, m_kd);
}

/******************************************************************/
/******************************************************************/
/**
 *  Reset distance, poles, and memory
 ******************************************************************/
void
CTransCoeff::reset(REAL rho, REAL kappa, int p)

{
	
  assert(p > 0 && p <= N_POLES);
	
  reallocate(p);
	
  m_p = p;
	
  initParams(rho, kappa);
  computeCoeff();
	
}

/******************************************************************/
/******************************************************************/
/**
 * Initialize the scalin factors for translation coefficient matrix.
 One corresponds to the radius of molecule and the other corresponds
 to the radius of the specific CG sphere w/in the molecule in the system
 ******************************************************************/
void
CTransCoeff::initScale(REAL scalei, REAL scaleo)
{
  m_scalei = scalei;
  m_scaleo = scaleo;
  m_fact.resize(N_POLES);
  m_fact[0] = 1.0;
  for(int n=1; n<N_POLES; n++)
		m_fact[n] = m_fact[n-1]*(m_scaleo / m_scalei);
}

/******************************************************************/
/******************************************************************/
/**
 *  Apply the translation operator to the MP coeffs
 ******************************************************************/
void
CTransCoeff::translate(const CMulExpan & Min, CLocalExpan & Mout,
											 int p, bool tpose) const
{
  assert(p <= N_POLES && Min.getRange().p2() >= p);
  Mout.reset(p);	
  if(tpose)
	{
		cout <<"not using transpose"<<endl; exit(1);
		
		assert( fabs (Min.getScale() - 1.0/m_scaleo) < 1e-5);
		Mout.setScale(m_scalei);		
		for (int n = 0; n < p; n++)
		{
			// m = 0 case : only real part
			int l;
			for (l = 0; l <= n; l++)
	      Mout(n,0) += TRANS(n,l,0)*Min(l,0); //real part
			for (; l < p; l++)
	    {
	      REAL T = (EVEN(n+l) ? TRANS(l,n,0) : -TRANS(l,n,0));
	      Mout(n,0) += FACT(l-n)*T * Min(l,0);
	    }      
			// for m > 0 part: have real and imaginary components
			for (int m = 1; m <= n; m++)
	    {
	      for (l = m; l <= n; l++)
				{
					Mout(n,2*m-1) += TRANS(n,l,m)*Min(l,2*m-1); //real part
					Mout(n,2*m  ) += TRANS(n,l,m)*Min(l,2*m  ); //imag part
				}
	      for (; l < p; l++)
				{
					REAL T = (EVEN(n+l) ? TRANS(l,n,m) : -TRANS(l,n,m));
					Mout(n,2*m-1) += FACT(l-n)*T * Min(l,2*m-1);
					Mout(n,2*m  ) += FACT(l-n)*T * Min(l,2*m  );
				}
	    }
		}
	}	
  // not transposed (ie forward)
  else
	{
		assert( fabs( Min.getScale() - 1.0/m_scalei ) < 1e-5);
		Mout.setScale(m_scaleo);
		
		for (int n = 0; n < p; n++)
		{
			// m = 0 case : only real part
			int l;
			for (l = 0; l <= n; l++)
	    {
	      REAL T = (EVEN(n+l) ? TRANS(n,l,0) : -TRANS(n,l,0));
	      Mout(n,0) += FACT(n-l)*T * Min(l,0);
	    }
			for (; l < p; l++)
	    {
	      Mout(n,0) += TRANS(l,n,0)*Min(l,0);
			}			
			// for m > 0 part: have real and imaginary components
			for (int m = 1; m <= n; m++)
	    {
	      for (l = m; l <= n; l++)
				{
					REAL T = (EVEN(n+l) ? TRANS(n,l,m) : -TRANS(n,l,m));
					Mout(n,2*m-1) += FACT(n-l)*T * Min(l,2*m-1);
					Mout(n,2*m  ) += FACT(n-l)*T * Min(l,2*m);
				}
	      for (; l < p; l++)
				{
					Mout(n,2*m-1) += TRANS(l,n,m)*Min(l,2*m-1);
					Mout(n,2*m  ) += TRANS(l,n,m)*Min(l,2*m  );
				}
	    }
		}    
	}	
}

/******************************************************************/
/******************************************************************/
/**
 *  Apply the derivative of the translation operator
 ******************************************************************/
void
CTransCoeff::dTranslate(const CMulExpan & Min, CLocalExpan & Mout,
												int p, bool tpose) const
{
  assert(p <= N_POLES && Min.getRange().p2() >= p);
  Mout.reset(p);
  if(tpose) Mout.setScale(m_scalei);
  else Mout.setScale(m_scaleo);
	
  for (int n = 0; n < p; n++)
	{
		// m = 0 case : only real part
		int l;
		if (tpose)
		{
			for (l = 0; l <= n; l++)
				Mout(n,0) += dTRANS(n,l,0)*Min(l,0); // real part
			
			for (; l < p; l++)
	    {
	      REAL dT = (EVEN(n+l) ? dTRANS(l,n,0) : -dTRANS(l,n,0));
	      Mout(n,0) += FACT(l-n)*dT * Min(l,0);
	    }
		}
		else
		{
			for (l = 0; l <= n; l++)
	    {
	      REAL dT = (EVEN(n+l) ? dTRANS(n,l,0) : -dTRANS(n,l,0));
	      Mout(n,0) += FACT(n-l)*dT * Min(l,0);
	    }
			
			for (; l < p; l++)
				Mout(n,0) += dTRANS(l,n,0)*Min(l,0);
		}
		
		// m > 0 part : have real and imaginary components
		for (int m = 1; m <= n; m++)
		{
			int l;
			if (tpose)
	    {
	      for (l = m; l <= n; l++)
				{
					Mout(n,2*m-1) += dTRANS(n,l,m)*Min(l,2*m-1); //real part
					Mout(n,2*m  ) += dTRANS(n,l,m)*Min(l,2*m  ); //imag part
				}
	      for (; l < p; l++)
				{
					REAL dT = (EVEN(n+l) ? dTRANS(l,n,m) : -dTRANS(l,n,m));
					Mout(n,2*m-1) += FACT(l-n)*dT * Min(l,2*m-1);
					Mout(n,2*m  ) += FACT(l-n)*dT * Min(l,2*m  );
				}
	    }
			else
	    {
	      for (l = m; l <= n; l++)
				{
					REAL dT = (EVEN(n+l) ? dTRANS(n,l,m) : -dTRANS(n,l,m));
					Mout(n,2*m-1) += FACT(n-l)*dT * Min(l,2*m-1);
					Mout(n,2*m  ) += FACT(n-l)*dT * Min(l,2*m);
				}
	      
	      for (; l < p; l++)
				{
					Mout(n,2*m-1) += dTRANS(l,n,m)*Min(l,2*m-1);
					Mout(n,2*m  ) += dTRANS(l,n,m)*Min(l,2*m  );
				}
	    }
		}
	}
}

/******************************************************************/
/******************************************************************/
/**
 * Function to compute translation coefficients
 ******************************************************************/
void
CTransCoeff::computeCoeff()
{
  TRANS(0,0,0) = m_exkid*m_K[0];
  m_exkid *= m_ir;
  
  if (m_bGrad)
    dTRANS(0,0,0) = -m_exkid*m_K[1];
  
  m_exkid *= m_scalei;  //was CMCoeff::RS;
	
  for (int l = 1; l < 2*m_p-1; l++)
	{
		TRANS(l,0,0) = m_exkid*m_K[l];
		m_exkid *= m_ir;
		
		if (m_bGrad)
			dTRANS(l,0,0) = m_exkid*(l*m_K[l] - (2*l+1)*m_K[l+1]);
		
		m_exkid *= m_scalei; //was CMCoeff::RS;
	}	
  computeCoeff_(m_T);  
  if (m_bGrad)
    computeCoeff_(m_dT);
}

/******************************************************************/
/******************************************************************/
/**
 *  Function to compute translation coefficients for l!=0.
 ******************************************************************/
void
CTransCoeff::computeCoeff_(vector<vector<vector<REAL> > > & U)
{
  // factors to rescale expansion from scalei to scaleo
  REAL kappasqr = CRExpan::KAPPA*CRExpan::KAPPA;
  REAL kpio = kappasqr * m_scalei * m_scaleo;
  REAL kpoo = kappasqr * m_scaleo * m_scaleo;
  REAL o_i  = m_scaleo / m_scalei;
  for (int l = 1; l < 2*m_p-2; l++)
	{
		U[l][1][0] = -(kpio*BETA(l-1,0)*U[l-1][0][0] +
									 o_i*ALPHA(l,0)*U[l+1][0][0])/ALPHA(0,0);
	}
  
  for (int n = 1; n < m_p-1; n++)
    for (int l = n+1; l < 2*m_p-2-n; l++)
		{
			U[l][n+1][0] = -(kpio*BETA(l-1,0)*U[l-1][n][0] +
											 kpoo*BETA(n-1,0)*U[l][n-1][0] +
											 o_i*ALPHA(l,0)*U[l+1][n][0])/ALPHA(n,0);
		}
	
  for (int m = 1; m < m_p; m++)
	{
		for (int l = m; l < 2*m_p-1-m; l++)
		{
			U[l][m][m] = -(kpio*DELTA(l,-m)*U[l-1][m-1][m-1] +
										 o_i*GAMMA(l+1,m-1)*U[l+1][m-1][m-1])/GAMMA(m,-m);
		}
		
		int n=m;
		for (int l = n+1; l < 2*m_p-2-n; l++)
		{
			U[l][n+1][m] = -(kpio*BETA(l-1,m)*U[l-1][n][m] +
											 o_i*ALPHA(l,m)*U[l+1][n][m])/ALPHA(n,m);
		}
		
		
		for (int n = m+1; n < m_p-1; n++)
			for (int l = n+1; l < 2*m_p-2-n; l++)
			{
				U[l][n+1][m] = -(kpio*BETA(l-1,m)*U[l-1][n][m] +
												 kpoo*BETA(n-1,m)*U[l][n-1][m] +
												 o_i*ALPHA(l,m)*U[l+1][n][m])/ALPHA(n,m);
			}
		
	}
}

/******************************************************************/
/******************************************************************/
/**
 * Computes the difference between
 inprod( tM2, S.tM1, p ) and inprod(tM2, S.tM1, p-1)
 ******************************************************************/
REAL
CTransCoeff::computeError(const CMulExpan & tM1, const CMulExpan & tM2, int p)
{
  assert(p <= m_p);
  assert(tM1.getRange().p2() >= p && tM2.getRange().p2() >= p);
	
  REAL err = 0.0;
	
  int l = p-1;
  for (int n = 0; n < p; n++)
	{
		err += TRANS(l,n,0)*tM1(n,0)*tM2(l,0);
		for (int m = 1; m <= n; m++)
			err += TRANS(l,n,m)*( tM1(n,2*m-1) * tM2(l,2*m-1) +
													 tM1(n,2*m)   * tM2(l,2*m));
	}
	
  int n = p-1;
  for (l = 0; l < p-1; l++)
	{
		REAL T = (EVEN(n+l) ? TRANS(n,l,0) : -TRANS(n,l,0));
		err += FACT(n-l)*T*tM1(n,0)*tM2(l,0);
	}
	
  for (int m = 1; m <= n; m++)
    for (l = m; l < p-1; l++)
		{
			REAL T = (EVEN(n+l) ? TRANS(n,l,m) : -TRANS(n,l,m));
			err += FACT(n-l)*T*( tM1(n,2*m-1)* tM2(l,2*m-1) +
													tM1(n,2*m)  * tM2(l,2*m));
		}
  
  return err;
}

/******************************************************************/
/******************************************************************/
/**
 *  Print out trans coefficients
 ******************************************************************/
void
CTransCoeff::outputTrans(int p) const
{
  cout << "****TRANS COEFFICIENTS****" << endl;
  for (int m = 0; m < p; m++)
	{
		cout << "\t---m = " << m << "---" << endl;
		for (int l = m; l < p; l++)
		{
			for (int n = m; n <= l; n++)
	    {
	      REAL r = fabs(TRANS(l,n,m))>1e-15 ? TRANS(l,n,m) : 0;
	      cout << r << " | ";
	    }
			cout << endl;
		}
	}
}

/******************************************************************/
/******************************************************************/
/**
 *  Print out derivatives of trans coefficients
 ******************************************************************/
void
CTransCoeff::outputdTrans(int p) const
{
  cout << "****dTRANS COEFFICIENTS****" << endl;
  for (int m = 0; m < p; m++)
	{
		cout << "\t---m = " << m << "---" << endl;
		for (int l = m; l < p; l++)
		{
			for (int n = m; n <= l; n++)
	    {
	      REAL r = fabs(dTRANS(l,n,m))>1e-15 ? dTRANS(l,n,m) : 0;
	      cout << r << " | ";
	    }
			cout << endl;
		}
	}
}

/******************************************************************/
/******************************************************************/
/**
 *  Export translation coefficient matrix as 2X2
 ******************************************************************/
void
CTransCoeff::exportMat(double ** T)
{
  int ind[m_p];
  ind[0] = 1;
  for (int j = 1; j < m_p; j++)
    ind[j] = ind[j-1] + 2*j-1;
	
  for (int n = 0; n < m_p; n++)
    for (int m = -n; m <= n; m++)
		{
			int j = ind[n] + (m+n);
			for (int l = abs(m); l < m_p; l++)
			{
				int k = ind[l] + (l+m);
				if (l > n)
					T[j][k] = TRANS(l,n,abs(m));
				else
					T[j][k] = (EVEN(n+l) ? TRANS(n,l,abs(m)) : -TRANS(n,l,abs(m)));
			}
		}
}

/******************************************************************/
/******************************************************************/
/**
 * Function to check the derivative rotation matrix size and
 ensure that there is a matrix coefficient at a given position.
 s_ = actual index address ( should be = n+s when called )
 ******************************************************************/
bool
CTransCoeff::checkT(int l, int n, int m)
{
  bool bN = (n >= 0 && n <= N_POLES-1 && n <=m_p-1);
  bool bL = (n <= l && l <= 2*m_p-2-n );
  bool bM = (m >= 0 && m <=n );
	
  if(!bN || !bL || !bM )
	{
		cout <<"Error with transcoeff T index access: "<<bN<<" "<<bL<<" "<<bM<<endl;
		cout <<"l="<<l<<" n="<<n<<" m="<<m<<"; m_p="<<m_p<<endl;
		return false;
	}  
  int sL = m_T.size();
  bool bsL=(l >= sL );
  if(bsL)
	{cout <<"Transcoeff T too small: l= "<<l<<" "<<sL<<endl; return false;}	
  int sN = m_T[l].size();
  bool bsN = (n >= sN );
  if(bsL)
	{cout <<"Transcoeff T too small: n="<<n<<" "<<sN<<endl; return false;}	
  int sM = m_T[l][n].size();
  bool bsM= (m >= sM );
  if(bsM)
	{cout <<"Transcoeff T too small: m="<<m<<" "<<sM<<endl; return false;}
	
  return true;
}
