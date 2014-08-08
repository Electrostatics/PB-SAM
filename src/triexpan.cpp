#include "triexpan.h"

/******************************************************************/
/******************************************************************/
/**
 * CExpCenter constructor.
 \param p an int of the number of poles of the triexpansion
 \param scale a scaling factor for the expansion trio
 ******************************************************************/
CTriExpan::CTriExpan(int p, REAL scale)
{
  assert(p >= 0);
  m_M[0].reset(p); m_M[1].reset(p); m_M[2].reset(p);
  m_M[0].setScale(scale); m_M[1].setScale(scale);m_M[2].setScale(scale);
}

/******************************************************************/
/******************************************************************/
/**
 * The printing operator for a triexpansion object
 ******************************************************************/
ostream &
operator<<(ostream & out, const CTriExpan & T)
{
  cout << "\t---000---" << endl;
  cout << T.m_M[0] << endl;
	
	
  cout << "\t---111---" << endl;
  cout << T.m_M[1] << endl;
	
  cout << "\t---222---" << endl;
  cout << T.m_M[2] << endl;
  
  return out;
}

/******************************************************************/
/******************************************************************/
/**
 * Output a complex number for a triexpansion class
 ******************************************************************/
void
CTriExpan::outputComplex(REAL fact) const
{
  for(int i=0; i< 3; i++)
    m_M[i].outputComplex(fact);
}

/******************************************************************/
/******************************************************************/
/**
 * Rotate a triexpansion object
 ******************************************************************/
void
CTriExpan::rotate(const CQuat & Q, int p)
{
  assert(p <= getRange().p2());
  CTriExpan G(p);
  G.setScale( getScale());
	
  for (int n = 0; n < p; n++)
	{
		// m = 0 case
		CPnt ar(m_M[0](n,0), m_M[1](n,0), m_M[2](n,0));
		CPnt br = Q*ar;
		G.m_M[0](n,0) = br.x();
		G.m_M[1](n,0) = br.y();
		G.m_M[2](n,0) = br.z();
		
		for (int m = 1; m <= n; m++)
	  {
	    CPnt ar(m_M[0](n,2*m-1), m_M[1](n,2*m-1), m_M[2](n,2*m-1));
	    CPnt ai(m_M[0](n,2*m),   m_M[1](n,2*m),   m_M[2](n,2*m));
	    CPnt br = Q*ar;
	    CPnt bi = Q*ai;
	    
	    G.m_M[0](n,2*m-1) = br.x();G.m_M[0](n,2*m) =  bi.x();
	    G.m_M[1](n,2*m-1) = br.y();G.m_M[1](n,2*m) =  bi.y();
	    G.m_M[2](n,2*m-1) = br.z();G.m_M[2](n,2*m) =  bi.z();
	  }
	}
  this->copy(G, p);
}

/******************************************************************/
/******************************************************************/
/**
 * rotate triexpan vectors, assume order has already been
 incremented elsewhere
 ******************************************************************/
void
CTriExpan::incRotate(const CQuat & Q)
{
  int p = getRange().p2();
  assert(p > 0);
  CTriExpan G(p, getScale());
	
  int n = p-1;
  // m = 0 case
  CPnt ar(m_M[0](n,0), m_M[1](n,0), m_M[2](n,0));
  CPnt br = Q*ar;
  G.m_M[0](n,0) = br.x();
  G.m_M[1](n,0) = br.y();
  G.m_M[2](n,0) = br.z();
	
  for (int m = 1; m <= n; m++)
	{
		CPnt ar(m_M[0](n,2*m-1), m_M[1](n,2*m-1), m_M[2](n,2*m-1));
		CPnt ai(m_M[0](n,2*m),   m_M[1](n,2*m),   m_M[2](n,2*m));
		CPnt br = Q*ar;
		CPnt bi = Q*ai;
		
		G.m_M[0](n,2*m-1) = br.x(); G.m_M[0](n,2*m) =  bi.x();
		G.m_M[1](n,2*m-1) = br.y(); G.m_M[1](n,2*m) =  bi.y();
		G.m_M[2](n,2*m-1) = br.z(); G.m_M[2](n,2*m) =  bi.z();
	}
  this->copy_p(G, p);
}


/******************************************************************/
/******************************************************************/
/**
 * CGradExpan class constructor.  Makes a gradient expansion class
 ******************************************************************/
CGradExpan::CGradExpan(const CPnt* g, /*const vector<CPnt> &g,*/ const vector<CPnt> &pos,
											 int p, bool bKappa, REAL scale, bool bMultipole)
{
  reset(p);
  setScale(scale);
	
  if(bMultipole)
    for(int i=0; i<pos.size(); i++)
		{
			if ( g[i].normsq() < 1e-15) continue;
			CSpPnt spos = CartToSph(pos[i]);
			CMExpan M(1.0, spos, bKappa, p, scale);
			m_M[0] += ( g[i].x() * M );
			m_M[1] += ( g[i].y() * M );
			m_M[2] += ( g[i].z() * M );
		}
  else
    for(int i=0; i<pos.size(); i++)
		{
			if ( g[i].normsq() < 1e-15) continue;
			CSpPnt spos = CartToSph(pos[i]);
			CLExpan L(1.0, spos, bKappa, p, scale);
			m_M[0] += ( g[i].x() * L );
			m_M[1] += ( g[i].y() * L );
			m_M[2] += ( g[i].z() * L );
		}
}
