#ifndef _EXPANSION_H_
#define _EXPANSION_H_

#include <vector>
#include <iostream>
#include <cassert>
#include "util.h"

/*###############################################################################
 * #
 * # File: expansion.h
 * #
 * # Date: June 2013
 * #
 * # Description: 
 * #
 * # Author: EH Yap, L Felberg
 * #
 * # Copyright ( c )
 * #
 * ################################################################################*/

#ifdef __DIFF__
	#define N_POLES 1
#else
	#define N_POLES 30
#endif

#ifdef __SSE
	#include "emmintrin.h"
#endif

#define PREC_LIMIT ( 1e-30 )

using namespace std;

/*#########################################################*/
/*#########################################################*/
////////////////////////////////////////////////
//    CRange
///////////////////////////////////////////////
/*#########################################################*/
/*#########################################################*/

class CRange
{
	// A range for an expansion [p1,p2 )
	// p1 is the index of the lowest order and p2 is one larger than the
	// index of the highest order.
	public:
		CRange(  ) : m_p1(0), m_p2(0) {}
		CRange( int p2 ) : m_p1(0), m_p2(p2) { makeValid(); }
		CRange( int p1, int p2 ) : m_p1(p1), m_p2(p2) { makeValid(); }

		// for calculating rotcoeffs
		CRange( int p1, int p2, bool bRot ) : m_p1(p1), m_p2(p2) 
			{ if( bRot ) assert(p2 <= 2*N_POLES-1); else makeValid(); } 

		const CRange & operator=( const CRange & R )
			{ m_p1 = R.p1(  ); m_p2 = R.p2(); return *this; }
		CRange( const CRange & R )
			{ *this = R; }

		bool operator==( const CRange & R ) const
			{ return ( (m_p1 == R.p1( )) && (m_p2 == R.p2())); }

		int p1(  ) const
			{ return m_p1; }
		int p2(  ) const
			{ return m_p2; }

		bool isEmpty(  ) const
			{ return ( m_p1 >= m_p2 ); }
		bool inRange( int n ) const
			{ return ( n >= m_p1 && n < m_p2 ); }
		void incRange(  )
			{ m_p2++; makeValid(  ); }
		void decRange(  )
			{ m_p2--; makeValid(  ); }

		// Compute the union of two ranges
		static CRange makeunion( const CRange & R1, const CRange & R2 )
		{
			int p1 = ( R1.p1( ) < R2.p1() ? R1.p1() : R2.p1());
			int p2 = ( R1.p2( ) > R2.p2() ? R1.p2() : R2.p2());
			return CRange( p1, p2 );
		}

		// Compute the intersection of two ranges
		static CRange intersection( const CRange & R1, const CRange & R2 ) 
		{
			int p1 = ( R1.p1( ) > R2.p1() ? R1.p1() : R2.p1());
			int p2 = ( R1.p2( ) < R2.p2() ? R1.p2() : R2.p2());
			return CRange( p1, p2 );
		}

		friend ostream & operator<<( ostream & out, const CRange & R );

	private:
		void makeValid(  )
		{
			if ( m_p1 >= m_p2 || m_p1 < 0 || m_p1 >= N_POLES )
			{
				m_p1 = 0;
				m_p2 = 0;
			}
			else if ( m_p2 > N_POLES )
				m_p2 = N_POLES;
		}
		int m_p1, m_p2;
};	// end class CRange

/*#########################################################*/
/*#########################################################*/
////////////////////////////////////////////////
//    CExpan
///////////////////////////////////////////////
/*#########################################################*/
/*#########################################################*/

class CExpan
{
	// A basic expansion.
	public:
		static void initConstants(  );

		// Choose and appropriate scale for the set of points
		static REAL chooseScale( const vector<CPnt> & pos );

		// Choose and appropriate scale given the minimum and maximum distances
		static REAL chooseScale( double min, double max );

		CExpan( int p = 0, REAL scale = 1.0 ) :  m_scale(scale) { reset(CRange(0,p)); }
		CExpan( const CRange & R, double scale = 1.0 ) : m_scale(scale)
			{ reset( R );}

		void copy( const CExpan & M, const int p );
		void copy_p( const CExpan & M, const int p );

		// Construct expansion from a vector of coefficients
		CExpan( const double * v, const CRange & R, double scale = 1.0 );
		CExpan( const vector<REAL> &M, const CRange & R, double scale = 1.0 );

		// Construct expansion from a row of a coefficients matrix
		// m - column length, i - index of row
		CExpan( const double * v, int m, int i, const CRange & R, double scale = 1.0 );

		// Copy constructor and assignment operator
		CExpan & operator=( const CExpan & M );
		CExpan( const CExpan & M )
			{ *this = M; }

		void recip(  )
			{ for ( int l = 0; l < m_M.size( ); l++) m_M[l] = -m_M[l]; }
		void reset( const CRange & R );
		void clear(  ) { for (int i = 0; i < length(); i++) m_M[i] = 0.0; }

		void setVector( const REAL *M, int len ) { m_M.assign(M, M+len);}

		// Accessing elements of expansion.
		REAL operator[]( int i ) const
			{ assert ( i < length( ) && i >= 0); return m_M[i]; }
		REAL & operator[]( int i ) 
			{ assert ( i < length( ) && i >= 0); return m_M[i]; }

		// Accessing the elements as complex numbers
		Complex comp( int n, int m ) const;

		// Accessing the elements by order and index
		REAL operator(  )(int n, int mm) const;
		REAL & operator(  )(int n, int mm);

		const vector<REAL> & getVector(  ) const
			{ return m_M; }

		const vector<REAL> getVector( int p ) const
		{ 
			vector<REAL> M; M.assign( m_M.begin( ), m_M.begin()+IDX[p]);
			return M;
		}

		// Accessing the scale of the expansion
		double getScale(  ) const
			{ return m_scale; }
		void setScale( double scale )
			{ m_scale = scale; }
		friend bool IsScaleEqual( const CExpan & M1, const CExpan & M2 );

		// Mathematical operators that change the expansion
		CExpan & operator*=( const REAL s );
		CExpan & operator*=( const REAL * C );
		CExpan & operator+=( const CExpan & E );
		CExpan & operator-=( const CExpan & E );
		CExpan operator-(  ) const;

		// Mathematical operations involving expansions
		friend CExpan operator+( const CExpan & E1, const CExpan & E2 );
		friend CExpan operator-( const CExpan & E1, const CExpan & E2 );
		friend double operator*( const CExpan & E1, const CExpan & E2 );
		friend CExpan operator*( const CExpan & E, const REAL * C )
			{ CExpan N = E; N *= C; return N; }
		friend CExpan operator*( REAL s, const CExpan & E );
		friend REAL inprod_unitScale( const CExpan & E1, const CExpan & E2 );

		// Accessing the range of the expansion.
		void setRange( const CRange & range ); 
		void setRange( int p )
			{ setRange( CRange(0, p )); }
		virtual const CRange & getRange(  ) const
			{ return m_range; }

		int length(  ) const
		{ return m_M.size(  ); }

		// Print-out functions
		friend ostream & operator<<( ostream & out, const CExpan & p );
		void output( int p1, int p2 );
		void output(  ) { output(m_range.p1(), m_range.p2()); }
		void outputComplex( REAL fact = 1.0 ) const;
		// Constants used in computing expansions
		static int IDX[2*N_POLES+1];
		static REAL SQRT2, ISQRT2;

		// Save and Undo functiona
		void saveUndo(  ) 
		{
			m_MU.assign( m_M.begin( ), m_M.begin()+IDX[m_range.p2()]
						-IDX[m_range.p1()] );
			m_rangeU = m_range;
			m_offsetU = m_offset;
			m_scaleU = m_scale;
		}
		void undo(  ) 
		{
			m_M.assign( m_MU.begin( ), m_MU.begin()+IDX[m_rangeU.p2()]
						-IDX[m_rangeU.p1()] );
			m_range = m_rangeU;
			m_offset = m_offsetU;
			m_scale = m_scaleU;
		}

		// debug
		bool isBlownup(  )
		{
			for ( int k=0; k<m_M.size( ); k++)
				if( fabs(m_M[k] ) > 5 || isnan(m_M[k])) return true;
				return false;
		}

		static REAL computeDev( const CExpan & M1, const CExpan & M2 );

	protected:
		static REAL RATIO;
		vector<REAL> m_M, m_MU;
		CRange m_range, m_rangeU;
		double m_scale, m_scaleU;

	private:
		int m_offset, m_offsetU;
};	// end class CExpan

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

class CLocalExpan;

/*#########################################################*/
/*#########################################################*/
////////////////////////////////////////////////
//    CMulExpan
///////////////////////////////////////////////
/*#########################################################*/
/*#########################################################*/

class CMulExpan : virtual public CExpan
{
	// A general multipole expansion class
	public:
		// Construct multipole expansion from a vector of coefficients
		CMulExpan( const double * v, const CRange & range, 
						double scale = 1.0 ) : 
		CExpan( v, range, 1.0/scale ) {}

		// Construct multipole expansion from another expansion's m_M
		// used now for converting a CExpan to a CMulExpan 
		CMulExpan( const vector<REAL> M, const CRange & range, double scale = 1.0 ) : 
		CExpan( M, range, scale ) {} // ***assume scale is already inversed

		// Construct the multipole expansion of a set of charges.
		CMulExpan( const vector<REAL> & ch, const vector<CPnt> & pos,
					int p, bool bKappa, REAL scale );
		/* CMulExpan( const vector<REAL> & ch, const vector<CPnt> & pos,
					int p, bool bKappa );*/

		// Create empty expansions
		CMulExpan( int p = 0 ) : CExpan(p) {}
		CMulExpan( const CRange & range, REAL scale = 1.0 ) 
						: CExpan(range, 1.0/scale) {}

		// Compute the inner product of a multipole and local expansions
		friend REAL inprod( const CMulExpan & M, const CLocalExpan & L );
};	// end class CMulExpan : virtual public CExpan

/*#########################################################*/
/*#########################################################*/
////////////////////////////////////////////////
//    CLocalExpan
///////////////////////////////////////////////
/*#########################################################*/
/*#########################################################*/

class CLocalExpan : virtual public CExpan
{
	// A general local expansion class
	public:
		// Construct the local expansion from a vector of coefficients
		CLocalExpan( const double * v, const CRange & range, 
						double scale = 1.0 ) : 
		CExpan( v, range, scale ) {}

		CLocalExpan( const vector<REAL> M, const CRange & range, 
						double scale = 1.0 ) : 
		CExpan( M, range, scale ) {}

		// Construct the local expansion of a set of charges
		CLocalExpan( const vector<REAL> & ch, const vector<CPnt> & pos,
						int p, bool bKappa, REAL scale );
		CLocalExpan( const REAL *ch, const vector<CPnt> & pos,
						int p, bool bKappa, REAL scale );

		/*
		CLocalExpan( const vector<REAL> & ch, const vector<CPnt> & pos,
		int p, bool bKappa );
		*/

		// Create empty expansions
		CLocalExpan( int p = 0 ) : CExpan(p) {}
		CLocalExpan( const CRange & range, REAL scale = 1.0 ) : CExpan(range,scale) {}

		// Compute the inner product of a multipole and local expansions
		friend REAL inprod( const CMulExpan & M, const CLocalExpan & L );
};	// end class CLocalExpan : virtual public CExpan

/*#########################################################*/
/*#########################################################*/
////////////////////////////////////////////////
//    CSHExpan
///////////////////////////////////////////////
/*#########################################################*/
/*#########################################################*/

class CSHExpan : virtual public CExpan
{
	// An expansion of spherical harmonics
	public:
		static void initConstants(  );

		CSHExpan(  ) : CExpan(), m_theta(0.0), m_phi(0.0) {}

		// Create the spherical harmonics up to order p for the given theta and phi.
		CSHExpan( REAL theta, REAL phi, int p ):
		CExpan( CRange(0,p )),  m_theta(theta),  m_phi(phi){computecoeff(p);}

		// Create a spherical harmonics expansion for calculating rotcoeffs
		// can go up to 2_NPOLES-1 and be incremented
		CSHExpan( REAL theta, REAL phi, int p, bool bRot ):
		CExpan( CRange(0,p,bRot )),  m_theta(theta),  m_phi(phi) 
			{ computecoeff( p, bRot ); }

		const CSHExpan & operator=( const CSHExpan & M )
		{ 
			m_theta = M.m_theta; 
			m_phi = M.m_phi;  
			m_legendre = M.m_legendre;
			static_cast<CExpan &>( *this ) = static_cast<const CExpan&>(M);
		}

		// Compute the derivative expansions
		CSHExpan dMdt(  ) const;
		CSHExpan dMdp(  ) const;

		// Increment and decrement the order of the expansion
		void inc(  );
		void dec(  ) {}

		// Compute a special degenerate expansion for the theta = 0.0 case
		static void specialSH( vector<REAL> & SH, int n, REAL val );
		CExpan m_legendre; // stores legendre for incrementing order later

	protected:
		// Convert the spherical harmonics to their derivative with repect to
		// either theta or phi.
		void derivTheta(  );
		void derivPhi(  );
		REAL m_theta, m_phi;

	private:
		// Computes the coefficients to the spherical harmonics expansion.
		void computecoeff( int p, bool bRot = false );

		// Compute the legendre polynomials up to order p
		void legendre(  );

		// Increment the order of the legendre polynomials
		void incLegendre(  );

		// Constants used to construct the spherical harmonics.
		static REAL CONST1[2*N_POLES][2*N_POLES], CONST2[2*N_POLES][2*N_POLES], 
					CONST3[2*N_POLES][2*N_POLES], CONST5[2*N_POLES],
					CONST6[2*N_POLES];
};	// end class CSHExpan : virtual public CExpan

/*#########################################################*/
/*#########################################################*/
////////////////////////////////////////////////
//    CRExpan
///////////////////////////////////////////////
/*#########################################################*/
/*#########################################################*/

class CRExpan : public CSHExpan
{
	// Radial expansion abstract base class ( adding the radial component of the
	// expansion to the spherical harmonics )
	public:
		static void initConstants( REAL kappa );

		const CRExpan & operator=( const CRExpan & M );

		// compute the derivative expansion with respect to rho.
		virtual CExpan dMdr(  ) const = 0;
		virtual void init(  )
			{  m_val = ( m_bKappa ? m_rho * KAPPA : 0.0 ); bessel(); }

		static REAL KAPPA;

	protected:
		// Construct a radial expansion for a single charge 
		CRExpan( REAL ch, const CSpPnt & spos, bool bKappa, int p ) :
				m_rho( spos.rho( )), m_bKappa(bKappa), m_k(ch), CExpan(p),
				CSHExpan( spos.theta( ), spos.phi(), p)     {}

		void applyRadialComponent(  );

		// Compute the appropriate bessel functions
		virtual void bessel(  )
			{  m_bessel.resize( m_range.p2( )+1, 1.0); }
		virtual void incBessel(  )
			{  m_bessel.resize( m_range.p2( )+1, 1.0); }

		static REAL CONST4[2*N_POLES];
		vector<REAL> m_bessel;
		bool m_bKappa;
		REAL m_rho, m_val, m_r, m_k;
};	// end class CRExpan : public CSHExpan

/*#########################################################*/
/*#########################################################*/
////////////////////////////////////////////////
//    CMExpan
///////////////////////////////////////////////
/*#########################################################*/
/*#########################################################*/

class CMExpan : public CRExpan, public CMulExpan
{
	// Multipole expansion for a single charge
	public:
		// Construct the multipole expansion of a single charge
		CMExpan( REAL ch,  const CSpPnt & spos, bool bKappa, int p,
				double scale = 1.0 ) : 
		CExpan( p, 1.0/scale ), CRExpan(ch, spos, bKappa, p)
			{ init(  ); applyRadialComponent();}
		CMExpan( const CMExpan & M ) : 
				CRExpan(static_cast<const CRExpan&>(M)) {}

		virtual void init(  );
		const CMExpan & operator=( const CMExpan & M )
		{ 
			// static_cast<CRExpan&>( *this ) 
			// 			= static_cast<const CRExpan&>(M);
			static_cast<CMExpan&>( *this ) = static_cast<const CMExpan&>(M);
			return *this;
		}

		// Compute the derivative with respect to the radius
		virtual CExpan dMdr(  ) const;

		static void BESSEL( vector<REAL> & K, REAL val );
		static void INCBESSEL( vector<REAL> & K, REAL val );

	protected:
		virtual void bessel(  );
		virtual void incBessel(  );  
};

/*#########################################################*/
/*#########################################################*/
////////////////////////////////////////////////
//    CLExpan
///////////////////////////////////////////////
/*#########################################################*/
/*#########################################################*/

class CLExpan : public CRExpan, public CLocalExpan
{
	// Local expansion for a single charge
	public:
		// Construct the local expansion of a single charge.
		CLExpan( REAL ch, const CSpPnt & spos, bool bKappa, int p, 
					double scale = 1.0 ) : CExpan( p, scale), 
					CRExpan(ch, spos, bKappa, p) 
			{ init(  ); applyRadialComponent(); }

		CLExpan( const CLExpan & L ) : CRExpan(static_cast<const CRExpan&>(L)) {}

		virtual void init(  );
		const CLExpan & operator=( const CLExpan & M )
		{ 
			static_cast<CRExpan&>( *this ) = static_cast<const CRExpan&>(M);
			return *this;
		}

		// Compute the derivative with respect to the radius
		virtual CExpan dMdr(  ) const;

		static void BESSEL( vector<REAL> & K, REAL val );
		static void INCBESSEL( vector<REAL> & K, REAL val );

	protected:
		virtual void bessel(  );
		virtual void incBessel(  );
};	// end class CLExpan : public CRExpan, public CLocalExpan

/*#########################################################*/
/*#########################################################*/
// A templated constructor for one 
// of the two types of expansions
/*#########################################################*/
/*#########################################################*/

template<class T> class CRExpanConstructor
{
	public:
	CRExpanConstructor( bool bKappa = false ) : m_bKappa(bKappa) {}
	T operator(  )(double ch, const CSpPnt & q, int p, double scale) const
		{ return T( ch, q, m_bKappa, p, scale ); }

	protected:
	bool m_bKappa;
}; // end class CRExpanConstructor

/*#########################################################*/
/*#########################################################*/
// Clear the expansion and set its range
/*#########################################################*/
/*#########################################################*/

inline void 
CExpan::reset( const CRange & R )
{
	m_M.resize( IDX[R.p2( )] - IDX[R.p1()]);
	m_M.assign( m_M.size( ), 0.0);

	// If the range does not start at 0 we save an offset value.
	m_offset = IDX[R.p1(  )];
	m_range = R;
}	// end CExpan::reset

/*#########################################################*/
/*#########################################################*/
// Convert the real valued expansion to complex form.
/*#########################################################*/
/*#########################################################*/

inline Complex
CExpan::comp( int n, int m ) const
{
	assert( n >= m_range.p1( ) && n < m_range.p2() && abs(m) <= n);

	if ( m == 0 )
		return Complex( (*this )(n,0), 0.0);
	else
	{
		int i = IDX[n] - m_offset + 2*abs( m )-1;
		if ( m > 0 )
			return ISQRT2*Complex( m_M[i], m_M[i+1] );
		else
			return ISQRT2*Complex( m_M[i], -m_M[i+1] );
	}
}	// end CExpan::comp( int n, int m ) const

/*#########################################################*/
/*#########################################################*/
// Return a coefficient of order n,mm ( 0 <= mm <= 2*n+1 )
/*#########################################################*/
/*#########################################################*/

inline REAL 
CExpan::operator(  )( int n, int mm ) const
{ 
	assert( n < m_range.p2() && mm < 2*n+1 ); 
	return m_M[IDX[n] /*- m_offset*/ + mm];
}	// end CExpan::operator(  )

/*#########################################################*/
/*#########################################################*/
// Return a reference to a coefficient 
// of order n,mm ( 0 <= mm <= 2*n+1 )
/*#########################################################*/
/*#########################################################*/

inline REAL & 
CExpan::operator(  )(int n, int mm) 
{ 
	assert( n < m_range.p2() && mm < 2*n+1); 
	return m_M[IDX[n] /*- m_offset*/ + mm];
}	// end CExpan::operator

/*#########################################################*/
/*#########################################################*/
// Set the range of an expansion. pad 
// with 0 if new coefficients are added.
/*#########################################################*/
/*#########################################################*/

inline void 
CExpan::setRange( const CRange & range )
{ 
	if ( range.p1( ) > m_range.p1() )
	{
		int minp2 = ( range.p2( ) < m_range.p2() 
						? range.p2() : m_range.p2());
		for ( int i = 0, j = IDX[range.p1( )]; j < minp2; j++, i++)
			m_M[i] = m_M[j];
		m_M.resize( IDX[range.p2( )] - IDX[range.p1()], 0.0);
	}
	else
	{
		m_M.resize( IDX[range.p2( )] - IDX[range.p1()], 0.0);
		int i = IDX[m_range.p1(  )]-IDX[range.p1()];
		for ( int j = 0; i < m_M.size( ), i != j; j++, i++ )
		{
			m_M[i] = m_M[j];
			m_M[j] = 0.0;
		}
	}
	m_range = range;
}	// end setRange

/*#########################################################*/
/*#########################################################*/
// Scale each level of the expansion by a 
// the corresponding element in C.
/*#########################################################*/
/*#########################################################*/

inline CExpan &
CExpan::operator*=( const REAL * C )
{
	vector<REAL>::iterator pT = m_M.begin(  );

	for ( int n = m_range.p1( ); n < m_range.p2(); n++ )
		for ( int m = 0; m < 2*n+1; m++ )
			*( pT++ ) *= C[n-m_range.p1()];

	return *this;
}	// end CExpan::operator*=

/*#########################################################*/
/*#########################################################*/
// Scale entire expansion by a scalar
/*#########################################################*/
/*#########################################################*/

inline CExpan &
CExpan::operator*=( const REAL s )
{
	vector<REAL>::iterator pT = m_M.begin(  );

	for ( int n = m_range.p1( ); n < m_range.p2(); n++)
		for ( int m = 0; m < 2*n+1; m++ )
			*( pT++ ) *= s;

	return *this;
}	// end CExpan::operator*=

/*#########################################################*/
/*#########################################################*/
// Add the coefficients of expansion M 
// to those of this expansion
/*#########################################################*/
/*#########################################################*/

inline CExpan &
CExpan::operator+=( const CExpan & M )
{
	//  if( !IsScaleEqual(*this, M )) {cout<<"scale error:" 
	//  <<this->getScale()<<" "<<M.getScale()<<" "<<M.getRange()<<endl;}
	assert(  IsScaleEqual(*this, M ) );

	CRange R = CRange::makeunion( m_range, M.m_range );
	setRange( R );

	int i = IDX[M.m_range.p1(  )] - IDX[m_range.p1()];
	int j = 0;
	for ( int j = 0; j < IDX[M.m_range.p2( )]
					-IDX[M.m_range.p1()]; j++, i++)
		m_M[i] += M.m_M[j];

	return *this;
}	// end CExpan::operator+=

/*#########################################################*/
/*#########################################################*/
// Subtract the coefficients of expansion 
// M from those of this expansion
/*#########################################################*/
/*#########################################################*/

inline CExpan &
CExpan::operator-=( const CExpan & M )
{
	assert(  IsScaleEqual(*this, M ) );
	CRange R = CRange::makeunion( m_range, M.m_range );
	setRange( R );

	int i = IDX[M.m_range.p1(  )] - IDX[m_range.p1()];
	for ( int j = 0; j < IDX[M.m_range.p2( )] 
				- IDX[M.m_range.p1()]; j++, i++)
		m_M[i] -= M.m_M[j];

	return *this;
}	// end CExpan::operator-=

/*#########################################################*/
/*#########################################################*/
// Make all coefficients equal to minus their value
/*#########################################################*/
/*#########################################################*/

inline CExpan
CExpan::operator-(  ) const
{
	CExpan E( *this );

	for ( int i = 0; i < length( ); i++ )
		E[i] = -E[i];

	return E;
}	// end CExpan::operator-

/*#########################################################*/
/*#########################################################*/
// Multiply the entire expansion by a scalar
/*#########################################################*/
/*#########################################################*/

inline CExpan 
operator*( REAL s, const CExpan & M )
{
	CExpan N( M );

	for ( int l = 0; l < M.length( ); l++ )
		N.m_M[l] *= s;

	return N;
} // end CExpan::operator*

/*#########################################################*/
/*#########################################################*/
// Compute the sum of two expansions
// The range of the result is the union of the two ranges
/*#########################################################*/
/*#########################################################*/

inline CExpan
operator+( const CExpan & M1, const CExpan & M2 )
{
	assert( IsScaleEqual(M1, M2 ));
	CRange R = CRange::makeunion( M1.getRange( ), M2.getRange());
	REAL scale = M1.getScale(  );
	CExpan S( R, scale );

	int i = CExpan::IDX[M1.m_range.p1(  )] - CExpan::IDX[R.p1()];
	for ( int j = 0; j < M1.length( ); j++, i++ )
		S.m_M[i] = M1.m_M[j];

	i = CExpan::IDX[M2.m_range.p1(  )] - CExpan::IDX[R.p1()];
	for ( int j = 0; j < M2.length( ); j++, i++ )
		S.m_M[i] = M2.m_M[j];

return S;
}	// end operator+

/*#########################################################*/
/*#########################################################*/
// Compute the difference of two expansions
// The range of the result is the union of the two ranges
/*#########################################################*/
/*#########################################################*/

inline CExpan
operator-( const CExpan & M1, const CExpan & M2 )
{
	assert( IsScaleEqual(M1, M2 ));
	CRange R = CRange::makeunion( M1.getRange( ), M2.getRange());
	REAL scale = M1.getScale(  );
	CExpan S( R, scale );

	int i = CExpan::IDX[M1.m_range.p1(  )] - CExpan::IDX[R.p1()];
	for  ( int j = 0; j < M1.length( ); j++, i++)
		S.m_M[i] = M1.m_M[j];

	i = CExpan::IDX[M2.m_range.p1(  )] - CExpan::IDX[R.p1()];
	for ( int  j = 0; j < M2.length( ); j++, i++)
		S.m_M[i] = M2.m_M[j];

	return S;
}	// end operator-

/*#########################################################*/
/*#########################################################*/
// Compute the inner product of a local 
// and a multipole expansion
/*#########################################################*/
/*#########################################################*/

inline REAL 
inprod( const CMulExpan & M,  const CLocalExpan & L )
{
	CRange R = CRange::intersection( M.getRange( ), L.getRange() );
	REAL sum( 0.0 );

	// In general: the appropriate scale is the product of the two scales
	// Our case: constructed such that lscale = 1/mscale -> fact = 1.0

	REAL fact = 1/( M.getScale( ) * L.getScale());

	if( fabs(fact-1.0 ) > 1e-5) 
		cout <<"fact "<< fact <<" "<<M.getScale(  )<<" "
					<< L.getScale()<<endl; 
	assert( fabs(fact-1.0 )< 1e-5 ); 

	REAL scale( 1.0 );
	for ( int n = 0; n < R.p1( ); n++)
		scale *= fact;

	for ( int n = R.p1( ); n < R.p2(); n++)
	{
		REAL s = 0.0;
		for ( int mm = 1; mm < 2*n+1; mm++ )
		{
			s += M( n,mm )*L(n,mm);
			//	cout<<n<<" "<<mm<<" "<<M( n,mm )<<" "<<L(n,mm)<<endl;
		}
		s *= 2.0; // accounts for negative m
		s += M( n,0 )*L(n,0); // m=0 case

		sum += scale * s;
		scale *= fact; 
	}
	return sum;
}	// end inprod

#ifdef __SSE
/*#########################################################*/
/*#########################################################*/
// Compute the inner product of two 
// expansions but ignore the scaling
// assume operation is performed on a unit sphere
/*#########################################################*/
/*#########################################################*/

	inline REAL 
	inprod_unitScale( const CExpan & E1,  const CExpan & E2 )
	{
		CRange R = CRange::intersection( E1.getRange( ), E2.getRange() );

		vector<REAL> ev1 = E1.getVector(  );
		vector<REAL> ev2 = E2.getVector(  );

		REAL sum( 0.0 );

		__m128d e1, e2;
		__m128d acc = _mm_setzero_pd(  );
		double temp[2];

		for ( int n = R.p1( ); n < R.p2(); n++)
		{
			int idxn = CExpan::IDX[n];
			for ( int mm = 1; mm < 2*n; mm+=2 ) 
			{
				e1 = _mm_loadu_pd(  &(ev1[idxn+mm] )  );
				e2 = _mm_loadu_pd(  &(ev2[idxn+mm] )  );
				acc = _mm_add_pd( acc, _mm_mul_pd(e1, e2 ));  
			}
		}

		_mm_storeu_pd( temp, acc );
		sum = 2.0* ( temp[0] + temp[1] );  // accounts for negative m

		for ( int n = R.p1( ); n < R.p2(); n++)
		{
			sum += E1( n,0 )*E2(n,0); // m=0 case
		}
		return sum;
	}	// end inprod_unitScale

	#else
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
	inline REAL 
	inprod_unitScale( const CExpan & E1,  const CExpan & E2 )
	{
		CRange R = CRange::intersection( E1.getRange( ), E2.getRange() );
		REAL sum( 0.0 );

		for ( int n = R.p1( ); n < R.p2(); n++)
		{
			REAL s = 0.0;
			for ( int mm = 1; mm < 2*n+1; mm++ ) s += E1(n,mm)*E2(n,mm);

			s *= 2.0; // accounts for negative m
			s += E1( n,0 )*E2(n,0); // m=0 case
			sum += s;
		}

		return sum;
	}	// end inprod_unitScale

#endif //endif-SSE

/*#########################################################*/
/*#########################################################*/
// Compute the inner product of a 
// local and a multipole expansion
/*#########################################################*/
/*#########################################################*/

inline double 
operator*( const CExpan & M1,  const CExpan & M2 )
{
	CRange R = CRange::intersection( M1.getRange( ), M2.getRange() );
	REAL sum( 0.0 );
	REAL fact = 1/( M1.getScale( ) * M2.getScale());
	/*
	if( fabs(fact-1.0 ) > 1e-5) 
	cout <<"fact "<< fact <<" "<<M.getScale(  )<<" "<< L.getScale()<<endl; 
	assert( fabs(fact-1.0 )< 1e-5 ); 
	*/

	REAL scale( 1.0 );
	for ( int n = 0; n < R.p1( ); n++)
		scale *= fact;

	for ( int n = R.p1( ); n < R.p2(); n++)
	{
		REAL s = 0.0;
		for ( int mm = 1; mm < 2*n+1; mm++ ){
			s += M1( n,mm )*M2(n,mm);
			//	cout<<n<<" "<<mm<<" "<<M( n,mm )<<" "<<L(n,mm)<<endl;
		}
		s *= 2.0; // accounts for negative m
		s += M1( n,0 )*M2(n,0); // m=0 case
		sum += scale * s;
		scale *= fact; 
	}
	return sum;
}	// end operator*

/*#########################################################*/
/*#########################################################*/
// corrected version of computeDev
/*#########################################################*/
/*#########################################################*/

inline REAL
CExpan::computeDev( const CExpan & M1, const CExpan & M2 )
{
	assert( M1.getRange( ) == M2.getRange());
	assert( IsScaleEqual(M1,M2 )); 

	REAL sum = 0;
	Complex tM1,tM2;

	int m_p = M1.getRange(  ).p2();

	for ( int n = 0; n < m_p; n++ )
	{
		for ( int m = 0; m <= n; m++ )
		{
			REAL s;

			tM1 = M1.comp( n,m );
			tM2 = M2.comp( n,m );

			if ( fabs(tM1.real( )) < PREC_LIMIT && 
						fabs( tM1.imag( )) < PREC_LIMIT)
			{	  
				if ( fabs(tM2.real( )) < PREC_LIMIT && 
							fabs( tM2.imag( )) < PREC_LIMIT)  
				{
					//cout <<"computedev1 "<<n<<" "<<m<<" "
					//			<<tM1<<" "<<tM2<<endl;
					continue;
				}
				else  
				{
					s = 1.0;
					//cout <<"computedev2 "<<n<<" "<<m<<" "
					//			<<tM1<<" "<<tM2<<endl;
				}
			}
			else 
			{
				if ( fabs(tM2.real( )) < PREC_LIMIT && 
							fabs( tM2.imag( )) < PREC_LIMIT)  
				{
					s = 1.0;
					//cout <<"computedev3 "<<n<<" "<<m<<" "
					//			<<tM1<<" "<<tM2<<endl;
				}
				else
				{
					//cout <<"computedev4 "<<n<<" "<<m<<" "
					// 			<<tM1<<" "<<tM2<<endl;
					Complex diff = ( tM1 - tM2 );
					REAL top = diff.real(  )*diff.real() +  diff.imag()*diff.imag();
					REAL bot = tM1.real(  )*tM1.real() +  tM1.imag()*tM1.imag() 
					+ tM2.real(  )*tM2.real() +  tM2.imag()*tM2.imag();

					if( m!=0 ) s = 2* top / bot; // for m!=0 
												//  to account for m<0 conjugates
					else s = top / bot;
					//cout <<"computedev4 "<<n<<" "<<m<<" "i
					//			<<tM1<<" "<<tM2<<" "<<s<<endl;
				}
			}
			sum += s;
		}	// end m
	}	// end n

	return sum;
}	// end computeDev

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

inline void 
CExpan::copy( const CExpan & M, const int p )
{
	assert( p <= M.m_range.p2( ));
	m_M.assign( M.m_M.begin( ), M.m_M.begin()+IDX[p]);
	m_range = CRange( 0,p,true ); // allows p up to 2NPOLES-1
	m_scale = M.m_scale;
}	// end copy

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

inline void 
CExpan::copy_p( const CExpan & M, const int p )
{
	assert(  p <= M.m_range.p2( ) && 
			( p == m_range.p2( ) || p == (m_range.p2()+1) ));
	m_M.resize( IDX[p-1] );
	m_M.insert( m_M.end( ), M.m_M.begin()+IDX[p-1], M.m_M.begin()+IDX[p]);
	m_range = CRange( 0,p, true ); // allows p up to 2NPOLES-1
	m_scale = M.m_scale;
}	// end copy_p

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

inline 
bool IsScaleEqual( const CExpan & M1, const CExpan & M2 )
	{ return(  fabs(M1.getScale( ) - M2.getScale()) < 1e-12);}
	//{ return (  M1.getScale( ) == M2.getScale() );}

#endif

/*#########################################################*/
/*#########################################################*/
// Commented out in old code!
/*#########################################################*/
/*#########################################################*/

/*
// Compute the inner product of a local and a multipole expansion
inline REAL 
inprod( const CMulExpan & M,  const CLocalExpan & L )
{
CRange R = CRange::intersection( M.getRange( ), L.getRange() );
REAL sum( 0.0 );

// The appropriate scale is the product of the two scales
REAL fact = 1/( M.getScale( ) * L.getScale());
if( fabs(fact-1.0 ) > 1e-5) 
cout <<"fact "<< fact <<" "<<M.getScale(  )<<" "<< L.getScale()<<endl; 
assert( fabs(fact-1.0 )< 1e-5 ); // fact should be 1.0

//cout <<"inprod: Scale_m = "<<M.getScale(  ) <<" Scale_l"  <<L.getScale()<<endl;
REAL scale( 1.0 );
for ( int n = 0; n < R.p1( ); n++)
scale *= fact;

for ( int n = R.p1( ); n < R.p2(); n++)
{
REAL s = 0.0;
for ( int mm = 0; mm < 2*n+1; mm++ ){
s += M( n,mm )*L(n,mm);
//	cout<<n<<" "<<mm<<" "<<M( n,mm )<<" "<<L(n,mm)<<endl;
}
sum += scale * s;
scale *= fact; 
}

return sum;
}
*/

/*
// wrong version ( Itay ) of computeDev : relative dev incorrectly 
// calculated

inline REAL
CExpan::computeDev( const CExpan & M1, const CExpan & M2 )
{
REAL sum = 0;
Complex tM1,tM2;
assert( M1.getRange( ) == M2.getRange());
assert( M1.getScale( ) == M2.getScale()); 
int m_p = M1.getRange(  ).p2();

for ( int n = 0; n < m_p; n++ )
for ( int m = 0; m <= n; m++ )
{
tM1 = M1.comp( n,m );
tM2 = M2.comp( n,m );
if ( tM1 == Complex(0.0,0.0 ) && tM2 == Complex(0.0,0.0))
continue;

Complex a;
if ( fabs(tM1.real( )) < PREC_LIMIT && 
fabs( tM1.imag( )) < PREC_LIMIT)
a = tM2;
else if ( fabs(tM2.real( )) < PREC_LIMIT && 
fabs( tM2.imag( )) < PREC_LIMIT)
a = tM1;
else
a = 0.5*( tM1 - tM2 )/(tM1 + tM2);

sum += a.real(  )*a.real() + a.imag()*a.imag();
}

return sum;
}
*/

