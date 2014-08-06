#ifndef _XFORMS_H_
#define _XFORMS_H_

#include <vector>
#include "util.h"
#include "expansion.h"
#include "expcenter.h"
#include "rotcoeff.h"
#include "transcoeff.h"
#include "triexpan.h"

#define MAX_REL_ERROR 1e-6
#define MIN_REL_ERROR ( MAX_REL_ERROR )
#define MAX_ERROR 1e-6

// For CGradExpan
#define dRHO 0
#define dTHETA 1
#define dPHI 2
#define dX 0
#define dY 1
#define dZ 2

class CMolecule;

//!  The CXFormBase class
/*!	The class that contains all details for an Abstract Base object
 for Reexpansion Transforms */
class CXFormBase
{
public:
	~CXFormBase(  ) {};
	
	static int computeOrder( REAL rho, REAL sourceTQ, REAL sourceScale );
	static void xformHrotMono( const CSolExpCenter &source, 
														const CSolExpCenter &target,CLocalExpan &Mout );
	static void xformHrotMonoG( const CSolExpCenter &source, 
														 const CSolExpCenter &target, CTriExpan &Mout );
	static void xformHrotMonoG( const CSolExpCenter &source, const CSolExpCenter &target, 
														 CTriExpan &Mout, double Q );
	static void xformGHrotMono( const CSolExpCenter &source, 
														 const CSolExpCenter &target, CTriExpan &Mout );
	static void xformGHrotMono( const CSolExpCenter &source, const CSolExpCenter &target, 
														 CTriExpan &Mout, CPnt gH );
	virtual void reset( const CPnt & P, int p ) {};
	
	virtual void xformF( const CMulExpan & Min, CLocalExpan & Mout, bool bFor ) {};
	virtual void xformF( const CTriExpan & Gin, CTriExpan & Gout, bool bFor ) {};
	virtual void xformH( const CMulExpan & Min, CLocalExpan & Mout, bool bFor ) {};
	virtual void xformH( const CExpan & Min, CGradExpan & Gout, bool bFor ) {};
	virtual void xformH( const CTriExpan & Gin, CTriExpan & Gout, bool bFor ) {};
	virtual int getOrder(  ) const {};
	
protected:
	CPnt m_P;
};

//!  The CXFormA class
/*!	The class that contains all details for an abstract object  */
class CXFormA : public CXFormBase
{
public:
	CXFormA(  ) :   m_pM1(NULL), m_pM2(NULL), m_rot(false) {}
	CXFormA( const CMulExpan & M1, const CMulExpan & M2, bool bGrad = false );
	CXFormA( const CSolExpCenter &C1, const CSolExpCenter &C2, bool bGrad = false );
	
	static void initConstants(  );
	
	virtual void reset( const CPnt & P, int p );
	virtual void xformH( const CMulExpan & Min, CLocalExpan & Mout, bool bFor );
	virtual void xformH( const CExpan & Min, CGradExpan & Gout, bool bFor );
	virtual void xformH( const CTriExpan & Gin, CTriExpan & Gout, bool bFor );
	virtual int getOrder(  ) const { return m_pH; }
	virtual void decOrderH(  ) {};
	
protected:
	CRotCoeff m_rot;
	CTransCoeff m_transH;
	int m_pH;
	
private:
	void sphToCart( CGradExpan & Gin );
	const CMulExpan *m_pM1, *m_pM2;
	CPnt m_R[3];
};

//!  The CXFormAIntra class
/*!	The class that contains all details for an intramolecular abstract object  */
class CXFormAIntra : public CXFormA
{
public: 
	void reset( const CPnt & P, int pH, int pF );
	CXFormAIntra( const CSolExpCenter & C1, const CSolExpCenter & C2 );

	virtual void xformF( const CMulExpan & Min, CLocalExpan & Mout, bool bFor );
	virtual void xformF( const CTriExpan & Gin, CTriExpan & Gout, bool bFor ) ;
	virtual void decOrderH(  ) {};

	void decOrderF(  );
	void incOrderF(  );
	int getOrderF(  ) const {return m_pF;}
	void setOrderF( int p );
	
private:
	CTransCoeff  m_transF;
	int m_pF;
};

//!  The CXFormN class
/*!	The class that contains all details for a numerical transform object  */
class CXFormN : public CXFormBase
{
public:
	CXFormN( const CSolExpCenter &C1, const CSolExpCenter &C2 );
	CXFormN( const CSolExpCenter &C1, const CMolecule &M2, REAL scale );
	
	virtual void reset( const CPnt & P, int p );
	virtual void xformH( const CMulExpan & Min, CLocalExpan & Mout, bool bFor );
	virtual void xformH( const CTriExpan & Gin, CTriExpan & Gout, bool bFor ) ;
	virtual int getOrder(  ) const {return N_POLES;} // HACK
	virtual double getScale1(  ) const {return m_scale1;}
	virtual double getScale2(  ) const {return m_scale2;}
	
protected:
	const CMulExpan *m_pM1, *m_pM2;
	int m_SPExSize1, m_SPExSize2;
	const CPnt *m_SPx1, *m_SPx2;
	const vector<double> &m_qSolvedH1, &m_qSolvedH2;
	const vector<CPnt> &m_gSolvedH1, &m_gSolvedH2;
	
	vector<CPnt> m_d1, m_d2;
	REAL m_relError;
	REAL m_scale1, m_scale2;
	const REAL *m_tQH1, *m_tQH2, *m_tGH1, *m_tGH2; 
	bool m_reset; //debug
	bool m_resetI, m_resetJ;
};

//!  The CXFormNIntra class
/*!	The class that contains all details for a intramolecular numerical transform object  */
class CXFormNIntra : public CXFormN
{
public:
	CXFormNIntra( const CSolExpCenter &C1, const CSolExpCenter &C2 );
	
	static void xformFH( const CSolExpCenter &source, const CSolExpCenter &target, 
											CLocalExpan & LFout, CLocalExpan & LHout, CPnt P, 
											const REAL* qSolvedF, const REAL* qSolvedH, 
											REAL tQF, REAL tQH );
	static void xformFH( const CSolExpCenter &source, const CSolExpCenter &target, 
											CLocalExpan & LFout, CLocalExpan & LHout, CPnt P );
	
	static void xformGFH( const CSolExpCenter &source, const CSolExpCenter &target, 
											 CTriExpan & LGFout, CTriExpan & LGHout, CPnt P );
	
	void xformF( const CMulExpan & Min, CLocalExpan & Mout, bool bFor );
	virtual void xformF( const CTriExpan & Gin, CTriExpan & Gout, bool bFor ) ;
	
private:
	const vector<double> &m_qSolvedF1, &m_qSolvedF2;
	const vector<CPnt>  &m_gSolvedF1, &m_gSolvedF2;
	const REAL *m_tQF1, *m_tQF2,*m_tGF1, *m_tGF2; 
};

#endif

