#include <cmath>
#include <cfloat>
#include "expansion.h"



#define EXPAN_EPS 1e-12

REAL CExpan::RATIO = 18.0; //pow(EXPAN_EPS, -2.0/N_POLES); // 4.0;
REAL CRExpan::KAPPA = 0.0;
int CExpan::IDX[2*N_POLES+1];


REAL CSHExpan::CONST1[2*N_POLES][2*N_POLES];
REAL CSHExpan::CONST2[2*N_POLES][2*N_POLES];
REAL CSHExpan::CONST3[2*N_POLES][2*N_POLES];
REAL CRExpan::CONST4[2*N_POLES];  
REAL CSHExpan::CONST5[2*N_POLES];
REAL CSHExpan::CONST6[2*N_POLES];

//REAL CExpan::SQRT2 = sqrt(2.0);
//REAL CExpan::ISQRT2 = sqrt(0.5);
REAL CExpan::SQRT2 = sqrt(1.0);
REAL CExpan::ISQRT2 = sqrt(1.0);


ostream & operator<<(ostream & out, const CRange & R)
{
  out << "[" << R.p1() << "," << R.p2() << "]" << endl;
  return out;
}

void
CRExpan::initConstants(REAL kappa)
{
  CSHExpan::initConstants();

  KAPPA = kappa;  

  for (int n = 0; n < 2*N_POLES; n++)
    CONST4[n] = 1.0/((2*n-1)*(2*n-3));
}

void
CExpan::initConstants()
{
  cout <<"N_POLES = "<<N_POLES<<endl;
  IDX[0] = 0;
  for (int n = 1; n <= 2*N_POLES; n++)
    IDX[n] = IDX[n-1]+2*n-1;
}

void
CSHExpan::initConstants()
{
  CExpan::initConstants();
  
  REAL temp[4*N_POLES];
  temp[0] = 1.0;
  for (int i = 1; i < 4*N_POLES; i++)
    temp[i] = temp[i-1]*sqrt((REAL)i);
  
  for (int n = 0; n < 2*N_POLES; n++)
    {
      for (int m = 0; m <= n; m++)
	{
	  CONST1[n][m] = (2*n-1)/(REAL)(n-m);
	  CONST2[n][m] = (n+m-1)/(REAL)(n-m);
	  CONST3[n][m] = temp[n-m]/temp[n+m];
	}

      CONST5[n] = (REAL)((2*n+1)*(2*n+3)); 
    }

  CONST6[0] = 1.0; CONST6[1] = 1.0;
  for (int n = 2; n < 2*N_POLES; n++)
    CONST6[n] = CONST6[n-1]*(2*n - 1);  
}

// Choose an appropriate scaling factor for the expansion based on min and max
// distance of the charges from the center
double
CExpan::chooseScale(double min, double max)
{
  double scale;
  if (max/min > RATIO)
    {
      cout << "BAD scale!!! min: " << min << " max: " 
	   << max << " ratio: " << max/min << endl;
      scale = min*sqrt(RATIO);
    }
  else
    scale = sqrt(max*min);

  return scale;
}

// Choose an appropriate scaling factor for the expansion based on the given
// set of charge positions
double
CExpan::chooseScale(const vector<CPnt> & pos)
{
  REAL mind = DBL_MAX, maxd = DBL_MIN;
  for (int i = 0; i < pos.size(); i++)
    {
      REAL d = pos[i].normsq();
      if (d < mind)
	mind = d;
      if (d > maxd)
	maxd = d;
    }
    
  return chooseScale(sqrt(mind), sqrt(maxd)); 
}

// Create an expansion using the supplied coefficients.
CExpan::CExpan(const double * v, const CRange & range, double scale) :
  m_scale(scale), m_offset(IDX[range.p1()])
{
  setRange(range);
  for (int i = 0; i < length(); i++)
    m_M[i] = v[i];
}

// Create an expansion using the supplied coefficients.
CExpan::CExpan(const vector<REAL> &M, const CRange & range, double scale) :
  m_scale(scale), m_offset(IDX[range.p1()])
{
  setRange(range);
  for (int i = 0; i < length(); i++)
    m_M[i] = M[i];
}

// Create an expansion using the supplied coefficients in the i-th row
// of the column-major matrix
CExpan::CExpan(const double * A, int m, int i, const CRange & range, 
	       double scale) : m_scale(scale), m_offset(IDX[range.p1()])
{
  setRange(range);
  for (int j = 0; j < length(); j++)
    m_M[j] = A[j*m+i];
}

CExpan & 
CExpan::operator=(const CExpan & M)
{  
  m_M = M.m_M; 
  m_range = M.m_range; 
  m_scale = M.m_scale;
  m_offset = M.m_offset;
  return *this;
}

// Construct a multipole expansion for a set of charges
CMulExpan::CMulExpan(const vector<REAL> & ch, const vector<CPnt> & pos,
		     int p, bool bKappa, REAL scale) 
  : CExpan(p)
{
  m_scale = 1.0/scale;

  for (int i = 0; i < pos.size(); i++)
    {
      if (ch[i] == 0.0)
	continue;

      CSpPnt spos = CartToSph(pos[i]);
      CMExpan M(ch[i], spos, bKappa, p, scale);
      *this += static_cast<CExpan&>(M);
    }


}

/*
// Construct a multipole expansion for a set of charges
CMulExpan::CMulExpan(const vector<REAL> & ch, const vector<CPnt> & pos,
		     int p, bool bKappa) : CExpan(p)
{
  double scale = CExpan::chooseScale(pos);
  for (int i = 0; i < pos.size(); i++)
    {
      if (ch[i] == 0.0)
	continue;

      CSpPnt spos = CartToSph(pos[i]);
      CMExpan M(ch[i], spos, bKappa, p, scale);
      *this += static_cast<CExpan&>(M);
    }

  m_scale = 1.0/scale;
}
*/

// Construct a local expansion for a set of charges
CLocalExpan::CLocalExpan(const REAL *ch, const vector<CPnt> & pos,
			 int p, bool bKappa, REAL scale) 
  : CExpan(p)
{
  m_scale = scale;

  for (int i = 0; i < pos.size(); i++)
    {
      if (ch[i] == 0.0)
	continue;

      CSpPnt spos = CartToSph(pos[i]);
      CLExpan L(ch[i], spos, bKappa, p, scale);
      *this += static_cast<CExpan&>(L);
    }

}

CLocalExpan::CLocalExpan(const vector<REAL> & ch, const vector<CPnt> & pos,
			 int p, bool bKappa, REAL scale) 
  : CExpan(p)
{
  m_scale = scale;

  for (int i = 0; i < pos.size(); i++)
    {
      if (ch[i] == 0.0)
	continue;

      CSpPnt spos = CartToSph(pos[i]);
      CLExpan L(ch[i], spos, bKappa, p, scale);
      *this += static_cast<CExpan&>(L);
    }

}
/*
// Construct a local expansion for a set of charges
CLocalExpan::CLocalExpan(const vector<REAL> & ch, const vector<CPnt> & pos,
			 int p, bool bKappa) : CExpan(p)
{
  double scale = CExpan::chooseScale(pos);
  for (int i = 0; i < pos.size(); i++)
    {
      if (ch[i] == 0.0)
	continue;

      CSpPnt spos = CartToSph(pos[i]);
      CLExpan L(ch[i], spos, bKappa, p, scale);
      *this += static_cast<CExpan&>(L);
    }

  m_scale = scale;
}
*/
const CRExpan & 
CRExpan::operator=(const CRExpan & M)
{
  m_bKappa = M.m_bKappa;
  m_rho = M.m_rho;
  m_val = M.m_val;
  m_r = M.m_r;
  m_k = M.m_k;
  m_bessel = M.m_bessel;
  static_cast<CSHExpan&>(*this) = static_cast<const CSHExpan&>(M); 
  return *this;
}

void 
CLExpan::init()
{
  CRExpan::init();
  
  m_r = 1.0/m_rho;
  m_k *= m_r*(m_bKappa ? exp(-m_val) : 1.0);
  m_r *= m_scale;

  //  cout <<"CLExpanInit "<<m_rho<<" "<<m_r<<" "<<m_bKappa<<" "<<m_val<<" "<<m_k <<endl;

}

// Apply the radial component of the expansion to the spherical harmonics
// part that was already computed.
void
CRExpan::applyRadialComponent()
{
  //  cout<<"m_r="<<m_r<<endl;
  //  cout<<"m_rho="<<m_rho<<endl;
  vector<REAL>::iterator it = m_M.begin();
  for (int n = 0; n < m_range.p2(); n++)
    {
      for (int m = 0; m < 2*n+1; m++){	
	*(it++) *= (m_bessel[n] * m_k);
      }
      //      cout << " radial: "<<n <<" "<<m_k<<endl;
      m_k *= m_r;
    }
}

void
CMExpan::init()
{
  CRExpan::init();
  m_r = m_rho*m_scale;
}

// Compute the Kirkwood MSBFs for the multipole expansion
void 
CMExpan::bessel() 
{
  CRExpan::bessel();
  CMExpan::BESSEL(m_bessel, m_val);
}

void
CMExpan::BESSEL(vector<REAL> & K, REAL val)
{
  if (val != 0.0)
    {
      REAL z = 0.5*val*val;
      for (int n = 0; n < K.size(); n++)
	{
	  REAL t = z/(2*n+3);
	  for (int j = 1; j <= 20; j++)
	    {
	      K[n] += t;
	      t *= (z/((j+1)*(2*(n+j)+3)));
	      if (t < 1e-20)
		break;
	    }
        }
    }
}

// Increment the order of the MSBFs by one
void 
CMExpan::incBessel() 
{
  CRExpan::incBessel();
  CMExpan::INCBESSEL(m_bessel, m_val);
}

void
CMExpan::INCBESSEL(vector<REAL> & K, REAL val)
{
  int p = K.size() + 1;
  K.resize(p);

  if (val != 0.0)
    {
      REAL z = 0.5*val*val;
      REAL t = z/(2*p + 3);
      for (int j = 1; j <= 20; j++)
	{
	  K[p] += t;
	  t *= (z/((j+1)*(2*(p+j)+3)));
	  if (t < 1e-20)
	    break;
	}
    }
}

// Compute the Kirkwood MSBFs for the local expansion
void 
CLExpan::bessel()
{
  CRExpan::bessel();
  CLExpan::BESSEL(m_bessel, m_val);
}

void
CLExpan::BESSEL(vector<REAL> & K, REAL val)
{
  if (val != 0.0)
    {
      if (K.size() > 1)
	{
	  K[1] += val; 

	  REAL valsqr = val*val;
	  for (int n = 2; n < K.size(); n++)
	    K[n] = K[n-1] + valsqr*K[n-2]*CONST4[n];
	}
    }

}

// Increment the order of the MSBFs by one
void 
CLExpan::incBessel()
{
  CRExpan::incBessel(); // increment range, pad with 1.0
  int p = m_range.p2();
  
  // calculate the p-1 coefficients
  if (p == 2)
    m_bessel[1] = m_val; 
  else
    {
      REAL valsqr = m_val*m_val;
      m_bessel[p-1] = m_bessel[p - 2] +
	valsqr*m_bessel[p - 3]*CONST4[p-1];
    }
}

void
CLExpan::INCBESSEL(vector<REAL> & K, REAL val)
{
  int p = K.size() + 1;
  K.resize(p);

  if (p == 2)
    K[p] = val; 
  else
    {
      REAL valsqr = val*val;
      K[p-1] = K[p - 2] + valsqr*K[p - 3]*CONST4[p-1];

    }
}

// Computes the coefficients to the spherical harmonics expansion.
void
CSHExpan::computecoeff(int p, bool bRot)
{

  // Generate the associated legendre polynomials
  legendre();

  // If this is used for rotcoeff, save the 
  // legendre in case we need to increment later
  if(bRot)
    {
      m_legendre.copy(*this, p);
      //     cout << "copied m_legendre "<<p<<" "<<m_legendre.getRange()<<endl;
    }
  
  REAL cosp[p], sinp[p], ang = 0;
  for (int m = 0; m < m_range.p2(); m++, ang += m_phi)
    {
      cosp[m] = cos(ang);
      sinp[m] = sin(ang);
    }
  
   vector<REAL>::iterator it = m_M.begin();
  for (int n = 0; n < m_range.p2(); n++)
    {
      *(it++) *= CONST3[n][0];
      REAL s = -SQRT2;
      for (int m = 1; m <= n; m++,s = -s)
	{
	  
	  *(it++) *= s*CONST3[n][m]*cosp[m];
	  *(it++) *= s*CONST3[n][m]*sinp[m];
	  
	}
    }
  
}


// Special spherical harmonics that are used in singular cases.
void 
CSHExpan::specialSH(vector<REAL> & SH, int n, REAL val)
{
  SH.resize(n);

  //  cout<<"SPECIAL SH ++++++++++++++++++++ "<<endl;
  //cout<<"n ="<<n<<" val = "<<val<<endl;
  if (n > 0)
    SH[0] = 0;
  if (n > 1)
    SH[1] = -1;
  if (n > 2)
    SH[2] = -val * 3;    

  for (int l = 3 ; l < n; l++)   
    {                             
      SH[l] = (val*(2*l-1)*SH[l-1] - l*SH[l-2])/(l-1);  
      //     cout<<"in a) " <<SH[l]<<endl;
    }
  // This part converts the legendre polynomial to a spherical harmonic.
  for (int l = 1; l < n; l++)
    {
      SH[l] *= (-CONST3[l][1]);
      // cout<<"in b) "<<CONST3[l][1]<<" "<<SH[l]<<endl;
    }
  /*  for (int l = 0; l < n; l++)
      cout<<SH[l]<<endl;*/

}

#ifdef __SSE
// Compute the legendre polynomials.
void 
CSHExpan::legendre() 
{                         
  REAL negterm = 1.0;
  REAL cost = cos(m_theta), sint = sin(m_theta);
  REAL sqrtterm = 1.0;   
  int p = m_range.p2();
  
  m_M[0] = CONST6[0];           
  negterm = -negterm;                                             
  sqrtterm *= sint;                                          
  if(p > 1)
    {                                                
      m_M[IDX[1]] = cost * m_M[0];           

     for(int n = 2; n < p; n++)
 
	m_M[IDX[n]] = cost*CONST1[n][0]*m_M[IDX[n-1]] - 
	  CONST2[n][0]*m_M[IDX[n-2]];}
  
  for (int m = 1, k = 1; m < p-1 ; m+=2, k+=4)
    {       

      int idxn = IDX[m];
      int idxn_plus1 = IDX[m+1];

      __m128d tmp1 = _mm_set_pd1( negterm*CONST6[m]*sqrtterm); 
      _mm_storeu_pd( &(m_M[idxn + k]), tmp1);
      __m128d tmp2 = _mm_set_pd1( -negterm*CONST6[m+1]*sqrtterm*sint ); 
      _mm_storeu_pd( &(m_M[idxn_plus1 + k+2]), tmp2);
      /*
      m_M[idxn       + k]   = negterm*CONST6[m]*sqrtterm;
      m_M[idxn       + k+1] = m_M[idxn + k];
      
      m_M[idxn_plus1 + k+2] = -negterm*CONST6[m+1]*sqrtterm*sint;
      m_M[idxn_plus1 + k+3] = m_M[idxn_plus1 + k+2];
      */

      sqrtterm *= sint*sint;                                          

      if (m < p-2)
	{                                                
	  int idxn_plus2 = IDX[m+2];
	
	  __m128d tmp1 = _mm_set_pd1( cost * (2*m+1) * m_M[idxn + k] ); 
	  _mm_storeu_pd( &(m_M[idxn_plus1 + k]), tmp1);
	  __m128d tmp2 = _mm_set_pd1( cost * (2*m+3) * m_M[idxn_plus1 + k+2] ); 
	  _mm_storeu_pd( &(m_M[idxn_plus2 + k+2]), tmp2);
	  __m128d tmp3 = _mm_set_pd1( cost*CONST1[m+2][m]*m_M[idxn_plus1 + k] -
				      CONST2[m+2][m]*m_M[idxn + k] ); 
	  _mm_storeu_pd( &(m_M[idxn_plus2 + k]), tmp3);
	  /*
	  m_M[idxn_plus1 + k]   = cost * (2*m+1) * m_M[idxn + k];  
	  m_M[idxn_plus1 + k+1] = m_M[idxn_plus1 + k];                    
	  m_M[idxn_plus2 + k+2] = cost * (2*m+3) * m_M[idxn_plus1 + k+2];  
	  m_M[idxn_plus2 + k+3] = m_M[idxn_plus2 + k+2];                    

	  m_M[idxn_plus2 + k] = cost*CONST1[m+2][m]*m_M[idxn_plus1 + k] -
	    CONST2[m+2][m]*m_M[idxn + k];
	  m_M[idxn_plus2 + k+1] =  m_M[idxn_plus2 + k];
	  */
	  for(int n = m+3; n < p; n++)
	    {
	      int idxn = IDX[n];
	      int idxn_minus1 = IDX[n-1];
	      int idxn_minus2 = IDX[n-2];

	      __m128d tmp1 = _mm_set_pd1( cost*CONST1[n][m]*m_M[idxn_minus1 + k] -
					  CONST2[n][m]*m_M[idxn_minus2 + k]  ); 
	      _mm_storeu_pd( &(m_M[idxn + k]), tmp1);
	      
	      __m128d tmp2 = _mm_set_pd1( cost*CONST1[n][m+1]*m_M[idxn_minus1 + k+2] -
					  CONST2[n][m+1]*m_M[idxn_minus2 + k+2] ); 
	      _mm_storeu_pd( &(m_M[idxn + k+2]), tmp2);
	      
	      /*
	      m_M[idxn + k] = cost*CONST1[n][m]*m_M[idxn_minus1 + k] -
                CONST2[n][m]*m_M[idxn_minus2 + k];
              m_M[idxn + k+1] =  m_M[idxn + k];

	      m_M[idxn + k+2] = cost*CONST1[n][m+1]*m_M[idxn_minus1 + k+2] -
                CONST2[n][m+1]*m_M[idxn_minus2 + k+2];
              m_M[idxn + k+3] =  m_M[idxn + k+2];
	      */
	    }                                                        
	}
                                                                
    }

  // compute omitted coeff
  if(p%2==0) // even
    {
      int n = p-1;
      int k = 2*n-1;
      __m128d tmp1 = _mm_set_pd1( negterm*CONST6[n]*sqrtterm  ); 
      _mm_storeu_pd( &(m_M[ IDX[n]+k ]), tmp1);
      /*
	m_M[IDX[n]+k]     = negterm*CONST6[n]*sqrtterm;
	m_M[IDX[n]+k+1]   = m_M[IDX[n]+k];
      */
    }
  else // odd
    {
      int n = p-2;
      int k = 2*n-1;

     __m128d tmp1 = _mm_set_pd1( cost * (2*n+1) * m_M[IDX[n]+k] ); 
      _mm_storeu_pd( &(m_M[ IDX[n+1]+k ]), tmp1);
      /*
      m_M[IDX[n+1]+k]   = cost * (2*n+1) * m_M[IDX[n]+k];  
      m_M[IDX[n+1]+k+1] = m_M[IDX[n+1]+k];       
      */             
    }

}

/*
// Compute the legendre polynomials.
void 
CSHExpan::legendre() 
{                         
  REAL negterm = 1.0;
  REAL cost = cos(m_theta), sint = sin(m_theta);
  REAL sqrtterm = 1.0;   
  int p = m_range.p2();
  
  m_M[0] = CONST6[0];           
  negterm = -negterm;                                             
  sqrtterm *= sint;                                          
  if(p > 1)
    {                                                
      m_M[IDX[1]] = cost * m_M[0];           

     for(int n = 2; n < p; n++)
 
	m_M[IDX[n]] = cost*CONST1[n][0]*m_M[IDX[n-1]] - 
	  CONST2[n][0]*m_M[IDX[n-2]];}
  
  for (int m = 1, k = 1; m < p-1 ; m+=2, k+=4)
    {       

      int n = m;
      m_M[IDX[n]+k]     = negterm*CONST6[n]*sqrtterm;
      m_M[IDX[n]+k+1]   = m_M[IDX[n]+k];
      m_M[IDX[n+1]+k+2] = -negterm*CONST6[n+1]*sqrtterm*sint;
      m_M[IDX[n+1]+k+3] = m_M[IDX[n+1]+k+2];
      sqrtterm *= sint*sint;                                          

      if (m < p-2)
	{                                                
	  m_M[IDX[n+1]+k]   = cost * (2*n+1) * m_M[IDX[n]+k];  
	  m_M[IDX[n+1]+k+1] = m_M[IDX[n+1]+k];                    
	  m_M[IDX[n+2]+k+2] = cost * (2*n+3) * m_M[IDX[n+1]+k+2];  
	  m_M[IDX[n+2]+k+3] = m_M[IDX[n+2]+k+2];                    

	  n=m+2;
	  m_M[IDX[n]+k] = cost*CONST1[n][m]*m_M[IDX[n-1]+k] -
	    CONST2[n][m]*m_M[IDX[n-2]+k];
	  m_M[IDX[n]+k+1] =  m_M[IDX[n]+k];
	  
	  for(int n = m+3; n < p; n++)
	    {

	      m_M[IDX[n]+k] = cost*CONST1[n][m]*m_M[IDX[n-1]+k] -
                CONST2[n][m]*m_M[IDX[n-2]+k];
              m_M[IDX[n]+k+1] =  m_M[IDX[n]+k];

	      m_M[IDX[n]+k+2] = cost*CONST1[n][m+1]*m_M[IDX[n-1]+k+2] -
                CONST2[n][m+1]*m_M[IDX[n-2]+k+2];
              m_M[IDX[n]+k+3] =  m_M[IDX[n]+k+2];

	    }                                                        
	}
                                                                
    }

  // compute omitted coeff
  if(p%2==0) // even
    {
      int n = p-1;
      int k = 2*n-1;
      m_M[IDX[n]+k]     = negterm*CONST6[n]*sqrtterm;
      m_M[IDX[n]+k+1]   = m_M[IDX[n]+k];
    }
  else // odd
    {
      int n = p-2;
      int k = 2*n-1;
      m_M[IDX[n+1]+k]   = cost * (2*n+1) * m_M[IDX[n]+k];  
      m_M[IDX[n+1]+k+1] = m_M[IDX[n+1]+k];                    
    }

}

*/
#else 
// Compute the legendre polynomials.
void 
CSHExpan::legendre() 
{                         
  REAL negterm = 1.0;
  REAL cost = cos(m_theta), sint = sin(m_theta);
  REAL sqrtterm = 1.0;   
  int p = m_range.p2();
  
  m_M[0] = CONST6[0];           
  negterm = -negterm;                                             
  sqrtterm *= sint;                                          
  if(p > 1)
    {                                                
      m_M[IDX[1]] = cost * m_M[0];           

     for(int n = 2; n < p; n++)
 
	m_M[IDX[n]] = cost*CONST1[n][0]*m_M[IDX[n-1]] - 
	  CONST2[n][0]*m_M[IDX[n-2]];}
  
  for (int m = 1, k = 1; m < p ; m++, k+=2)
    {       
      m_M[IDX[m]+k] = negterm*CONST6[m]*sqrtterm;
      m_M[IDX[m]+k+1] =  m_M[IDX[m]+k];
      negterm = -negterm;                                             
      sqrtterm *= sint;                                          

      if (m < p-1)
	{                                                
	  m_M[IDX[m+1]+k] = cost * (2*m+1) * m_M[IDX[m]+k];  
	  m_M[IDX[m+1]+k+1] =  m_M[IDX[m+1]+k];                    

	  for(int n = m+2; n < p; n++)
	    {
	      m_M[IDX[n]+k] = cost*CONST1[n][m]*m_M[IDX[n-1]+k] -
                CONST2[n][m]*m_M[IDX[n-2]+k];
              m_M[IDX[n]+k+1] =  m_M[IDX[n]+k];

	    }                                                        
	}                                                                
    }
}

#endif // endif-SSE

// calculate the legendre polynomials for the new p-l layer (to m_legendre)
// then copy the new layer to m_M
void
CSHExpan::incLegendre()
{
  int p = m_range.p2();
  int n = p - 1;
  int m;
  REAL cost = cos(m_theta), sint = sin(m_theta);

  assert(n >= 1);
  
  m_legendre.setRange( CRange(0,p,true) );

  if(n == 1)
    {
      m = 1;  
      REAL negterm = -1.0;
      REAL sqrtterm = sint;   
      m_legendre[IDX[m]+2*m-1] = negterm*CONST6[m]*sqrtterm;    
      m_legendre[IDX[m]+2*m] =  m_legendre[IDX[m]+m];

      m = 0;
      m_legendre[IDX[m+1]] = cost * (2*m+1) * m_legendre[IDX[m]+m];  
    }

  else
    {
      
      // recursive propagation from previous n's need to 
      // use m_legendre since m_Ms are modified in computecoeff

      m_legendre[IDX[n]] = cost*CONST1[n][0]*m_legendre[IDX[n-1]] - 
	CONST2[n][0]*m_legendre[IDX[n-2]];
  
      m = n;  
      REAL negterm = ( n % 2 == 0? 1.0 : -1.0);
      REAL sqrtterm = pow( sint, m); 
      m_legendre[IDX[m]+2*m-1] = negterm*CONST6[m]*sqrtterm;    
      m_legendre[IDX[m]+2*m] =  m_legendre[IDX[m]+2*m-1];
      
      m = n-1;
      m_legendre[IDX[m+1]+2*m-1] = cost * (2*m+1) * m_legendre[IDX[m]+2*m-1];  
      m_legendre[IDX[m+1]+2*m]   = m_legendre[IDX[m+1]+2*m-1];                    
      
      for (int m = 1; m < n-1 ; m++)
	{       
	  m_legendre[IDX[n]+2*m-1] = cost*CONST1[n][m]*m_legendre[IDX[n-1]+2*m-1] - 
	    CONST2[n][m]*m_legendre[IDX[n-2]+2*m-1];
	  m_legendre[IDX[n]+2*m] =  m_legendre[IDX[n]+2*m-1];  
	}
    }

  // update m_legendre for future increments
      (*this).copy_p(m_legendre, p);
      //cout <<"copyed leg to M: p "<<p<<" range M "<<m_range<<endl;
      
}


// Generate the derivative of the expansion with respect to the radius
CExpan
CMExpan::dMdr() const
{
  double E[length()];
  double * pE = &(E[0]);
  double ir = 1.0/m_rho;

  vector<REAL>::const_iterator it = m_M.begin();
  REAL c = (m_bKappa ? KAPPA*KAPPA*m_rho*m_rho : 0.0);
  for (int n = m_range.p1(); n < m_range.p2(); n++)
    {
      REAL C = ir*(n + (m_bKappa ? c*m_bessel[n+1]/(m_bessel[n]*(2*n+3)) 
			: 0.0));
      for (int m = 0; m < 2*n+1; m++)
	*(pE++) = *(it++) * C;
    }

  return CExpan(E, m_range, getScale());
}

// Generate the derivative of the expansion with respect to the radius
CExpan
CLExpan::dMdr() const
{
  double E[length()];
  double * pE = &(E[0]);
  double ir = 1.0/m_rho;

  vector<REAL>::const_iterator it = m_M.begin();
  for (int n = m_range.p1(); n < m_range.p2(); n++)
    {
      REAL C = ir*(n - (m_bKappa ? (2*n+1)*m_bessel[n+1]/m_bessel[n] : 0.0));

      for (int m = 0; m < 2*n+1; m++)
	*(pE++) = *(it++) * C;
    }

  return CExpan(E, m_range, getScale());
}

void 
CSHExpan::inc()
{
  int p = m_range.p2() + 1; // new no. of poles
  assert(p <= 2*N_POLES);
  setRange( CRange(0, p, true) ); // note that we need this format of CRange circumvent makeValid()

  incLegendre();
  //cout <<" after legendre "<<endl;
  //  this->outputComplex(); 

  REAL cosp[p], sinp[p], ang = 0;
  for (int m = 0; m < p; m++, ang += m_phi)
    {
      cosp[m] = cos(ang);
      sinp[m] = sin(ang);
    }

  m_M[IDX[p-1]] *=  CONST3[p-1][0];

  REAL s = -SQRT2;
  for (int m = 1; m <= p-1; m++, s = -s)
    {
      m_M[IDX[p-1] + 2*m-1]   *=  s*CONST3[p-1][m]*cosp[m];
      m_M[IDX[p-1] + 2*m]     *=  s*CONST3[p-1][m]*sinp[m];
    }

}





CSHExpan
CSHExpan::dMdt() const
{
  CSHExpan M(*this);
  M.derivTheta();

  return M;
}

// Convert the spherical harmonics expansion to its derivative with respect
// to theta.
void
CSHExpan::derivTheta()
{
  REAL cot_t = 1.0/tan(m_theta);
  REAL sinp = sin(m_phi), cosp = cos(m_phi);

  int j = 0;
  vector<REAL>::iterator pT = m_M.begin();
  
  *(pT++) = 0.0; j++;
  for (int n = 1; n < m_range.p2(); n++)
    {
      j++;
      REAL k1 = sqrt((REAL)n*(n+1));
      *(pT++) = -k1 * ISQRT2*(cosp*m_M[j] + sinp*m_M[j+1]);
     
      REAL k2 = cot_t;
      for (int m = 1; m < n; m++)
	{
	  k1 = sqrt((REAL)(n-m)*(n+m+1));
	  *(pT++) = k2*m_M[j] - k1*(cosp*m_M[j+2] + sinp*m_M[j+3]);

	  *(pT++) = k2*m_M[j+1] - k1*(cosp*m_M[j+3] - sinp*m_M[j+2]);

	  k2 += cot_t;
	  j += 2;
	}

      *(pT++) = k2*m_M[j]; *(pT++) = k2*m_M[j+1];
      j += 2;
    }
}

CSHExpan
CSHExpan::dMdp() const
{
  CSHExpan M(*this);
  M.derivPhi();

  return M;
}

// Convert the spherical harmonics expansion to its derivative with respect
// to phi.
void
CSHExpan::derivPhi() 
{
  vector<REAL>::const_iterator pS = m_M.begin();
  vector<REAL>::iterator pT = m_M.begin();

  for (int n = 0; n < m_range.p2(); n++)
    {
      *(pT++) = 0.0;
     
      for (int m = 1; m <= n; m++)
      {
	vector<REAL>::iterator pS = pT;
	pS++;
	REAL tmp = *pT;
	*(pT++) = -m*(*pS);
	*(pT++) = m*tmp;
      }
    }
}

ostream & 
operator<<(ostream & fout, const CExpan & M)
{
  vector<REAL>::const_iterator pS = M.m_M.begin();

  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  for (int n = M.m_range.p1(); n < M.m_range.p2(); n++)
    {
      for (int m = 0; m < 2*n+1; m++) 
	{
	  REAL out = *(pS++);
	  //if(abs(out) < 1e-5 ) out = 0.0;
	  if(abs(out) < 1e-15 ) out = 0.0;
	  fout << out << " ";
	}
      fout << endl;
    }

  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  return fout;
}

void
CExpan::output(int p1, int p2)
{
  CRange temp = m_range;
  m_range = CRange(p1,p2);
  cout << *this; 
  m_range = temp;
}

void
CExpan::outputComplex(REAL fact) const
{
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  for (int n = m_range.p1(); n < m_range.p2(); n++)
    {
      for (int m = 0; m <=n; m++) 
	{
	  Complex out = comp(n,m);
	  REAL re =( (abs(out.real()) < 1e-15 )? 0.0 : out.real());
	  REAL im =( (abs(out.imag()) < 1e-15 )? 0.0 : out.imag());
	  cout << "("<<fact*re<< ","<<fact*im<<") " ;
	}
      cout << endl;
    }

  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

}


