#ifndef _UTIL_H_
#define _UTIL_H_

#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#include <memory>

using namespace std;

typedef double REAL;
typedef complex<REAL> Complex;

class CPnt;
class CQuat;

//!  The CSpPnt class 
/*! The spherical coordinate class. Contains a set of
spherical coordinates and operations that may be performed on it. */
class CSpPnt
{
	public:
	//!  The CSpPnt class constructor
	/*! The spherical coordinate class constructor. Initializes
	 a spherical coordinate class of (0,0,0)  */
	CSpPnt(  ) 
	{ zero(  ); }
	//!  The CSpPnt class constructor
	/*! The spherical coordinate class constructor: (r, t, p). 
	 \param rho an input r for radius
	 \param theta an input t angle
	 \param phi an input p angle */
	CSpPnt( REAL rho_, REAL theta_, REAL phi_ ) :
	m_rho( rho_ ), m_theta(theta_), m_phi(phi_) {}
	
	//!  The CSpPnt class = operator
	const CSpPnt & operator=( const CSpPnt & s ) 
	{ m_rho = s.rho(  ); m_theta = s.theta(); m_phi = s.phi(); return *this; }
	
	//!  The CSpPnt class zero function
	void zero(  )
	{ m_rho = 0.0; m_theta = 0.0; m_phi = 0.0; }
	
	// Printing functions
	//!  The CSpPnt class print rho
	REAL & rho(  ) 
	{ return m_rho; }
	//!  The CSpPnt class print theta
	REAL & theta(  ) 
	{ return m_theta; }
	//!  The CSpPnt class print phi
	REAL & phi(  ) 
	{return m_phi; }
	//!  The CSpPnt class print rho
	const REAL rho(  ) const
	{ return m_rho; }
	//!  The CSpPnt class print theta
	const REAL theta(  ) const
	{ return m_theta; }
	//!  The CSpPnt class print phi
	const REAL phi(  ) const
	{return m_phi; }
	
	//!  The CSpPnt class << operator to print out rho, theta, phi
	friend ostream & operator<<( ostream & out, const CSpPnt & p );
	
	//!  The CSpPnt class convert spherical to cartesian
	/*!  \param s a spherical coordinate object to convert
			 \return a cartesian coordinate object */
	friend CPnt SphToCart( const CSpPnt & s );

	//!  The CSpPnt class convert cartesian to spherical
	/*!  \param c a cartesian coordinate object to convert
			 \return a spherical coordinate object */
	friend CSpPnt CartToSph( const CPnt & c );
	
private:
	REAL m_rho;				//!< R part of a spherical coordinate
	REAL m_theta;			//!< Theta part of a spherical coordinate
	REAL m_phi;				//!< Phi part of a spherical coordinate
};

REAL dot( const CPnt & p1, const CPnt & p2 );

//!  The CPnt class 
/*! The cartesian coordinate class. Contains a set of
cartesian coordinates and operations that may be performed on it. */
class CPnt
{
public:
	//!  The CPnt class constructor
	/*! The cartesian coordinate class constructor. Initializes
	 a cartesian coordinate class of (0,0,0)  */
	CPnt(  ) 
	{ zero(  ); }
		//!  The CPnt class constructor
	/*! The cartesian coordinate class constructor: (x,y,z). 
	 \param x an input x
	 \param y an input y
	 \param z an input z */
	CPnt( REAL x_, REAL y_, REAL z_ ) :
	m_x( x_ ), m_y(y_), m_z(z_) {}
	//!  The CPnt class constructor
	/*! The cartesian coordinate class constructor: (x,y,z). 
	 \param p an input cartesian coord object */
	CPnt( const CPnt & p ) : 
	m_x( p.x( )), m_y(p.y()), m_z(p.z()) {}
	
	//!  The CPnt class = operator	
	const CPnt & operator=( const CPnt & c )
	{ m_x = c.x(  ); m_y = c.y(); m_z = c.z(); return *this; }

	//!  The CPnt class zero function
	void zero(  )
	{ m_x = 0.0; m_y = 0.0; m_z = 0.0; }

	//!  The CSpPnt class normsq function	
	/*!  Return the norm^2 of a vector ( i^2 + j^2 + k^2 )  */
	REAL normsq(  ) const
	{ return dot( *this,*this ); }
	//! The CPnt norm function
	/*!  Return the norm of a vector ( i^2 + j^2 + k^2 )^0.5  */
	REAL norm(  ) const
	{ return sqrt( normsq( )); }
	//! The CPnt normalize function
	/*!  Return a normalized vector ( vec / || vec || )  */
	const CPnt & normalize(  )
	{ ( *this ) *= (1.0/norm()); return *this; }
	//! The CPnt normalize function
	/*!  Return a normalized vector ( vec / || vec || )  */
	CPnt normalize(  ) const
	{ CPnt P = *this; P *= ( 1.0/norm( )); return P; }
	
// Inline functions, simple operators
	void operator+=( const CPnt & p )
	{ m_x += p.x(  ); m_y += p.y(); m_z += p.z(); }
	void operator-=( const CPnt & p )
	{ m_x -= p.x(  ); m_y -= p.y(); m_z -= p.z(); }
	void operator*=( REAL s )
	{ m_x *= s; m_y *= s; m_z *= s; }
	void operator/=( REAL s )
	{ REAL t = 1.0/s; *this *= t; }

	// Printing functions
	//!  The CSpPnt class print x
	const REAL x(  ) const
	{ return m_x; }
	//!  The CSpPnt class print y
	const REAL y(  ) const
	{ return m_y; }
	//!  The CSpPnt class print z
	const REAL z(  ) const
	{ return m_z; }
	//!  The CSpPnt class print x
	REAL & x(  ) 
	{ return m_x; }
	//!  The CSpPnt class print y
	REAL & y(  ) 
	{ return m_y; }
	//!  The CSpPnt class print z
	REAL & z(  ) 
	{ return m_z; }
	
// Inline functions described below
	friend CPnt operator+( const CPnt & p1, const CPnt & p2 );
	friend CPnt operator-( const CPnt & p1, const CPnt & p2 );
	
	friend CPnt cross( const CPnt & p1, const CPnt & p2 );
	friend CPnt operator*( REAL c, const CPnt & p1 );
	friend CPnt operator-( const CPnt & p );
	
	friend ostream & operator<<( ostream & out, const CPnt & p );
	friend CPnt SphToCart( const CSpPnt & s );
	friend CSpPnt CartToSph( const CPnt & c );
	friend REAL torsion( const CPnt & p1, const CPnt & p2,
											const CPnt & p3, const CPnt & p4 );
	friend REAL torsion( const CPnt & v1, const CPnt & v2,
											const CPnt & v3 );
	
private:
	REAL m_x;  //!< The x coordinate of the CPnt object
	REAL m_y;  //!< The y coordinate of the CPnt object
	REAL m_z;  //!< The z coordinate of the CPnt object
};

//!  The quaternion class
/*!
		The class that contains all details for a quaternion object. 
		Quaterions can be used to describe spatial rotations in space.
		It is represented here as:
		
		Q = m_real + m_imag.x*i + m_imag.y*j + m_imag.z*k
		
		the conjugate is 
		
		Q.conj = m_real - m_imag.x*i - m_imag.y*j - m_imag.z*k
		
		For more, see http://mathworld.wolfram.com/Quaternion.html
*/
class CQuat
{
public:
	 //!  The quaternion constructor
	/*!  Create an identity quaternion object  */	
	CQuat(  )
	{ identity(  ); }
	
	//!  The quaternion constructor
	/*!  Create an identity quaternion object from real and imag components
	 and normalize them.
	 \param real the real component of the quaternion
	 \param imag a CPnt object of imaginary components for the quaternion */	
	CQuat( REAL real, const CPnt & imag ) :
	m_real( real ), m_imag(imag) { normalize(); }
	
	//!  The quaternion constructor
	/*!  Create an identity quaternion object from an axis and an angle
	 \param axis a vector defining an axis of rotation
	 \param angle an angle to rotate around the axis by in radians 	 */	
	CQuat( const CPnt & axis, REAL angle )
	{  
		if( angle != 0.0 && axis.normsq( ) != 0.0)
		{
			m_real = cos( 0.5*angle );
			m_imag = sin( 0.5*angle )*(axis.normalize());
			normalize(  ); 
		}
		else identity(  );
	}
	
	//!  The quaternion constructor
	/*!  Create an identity quaternion object from another quaternion
	 \param q an input quaternion to create the new object from	 */		
	CQuat( const CQuat & q )
	{ *this = q; }
	
	CQuat & operator=( const CQuat & q )
	{ m_imag = q.m_imag; m_real = q.m_real; return *this; }
	
	//! CQuat identity
	/*! Generate the identity quaternion, with imag parts set to (0,0,0)
	 and the real part to 1.0 */	
	void identity(  )
	{ m_imag.zero(  ); m_real = 1.0;}
	
	//! CQuat normalize
	/*!  Normalize the quaternion, compute the norm and then
	 reset the imaginary and real parts of the quaternion  */
	void normalize(  )
	{ 
		REAL in_norm = 1.0/sqrt( m_imag.normsq( ) + m_real*m_real);
		m_imag *= in_norm; m_real *= in_norm;
	}  
	//! CQuat conj
	/*! Set the imaginary part to their conjugates. */
	void conj(  )
	{ m_imag = -m_imag; }
	
	// Printing the different parts of the quaternion
	REAL real(  ) const
	{ return m_real; }
	const CPnt & imag(  ) const
	{ return m_imag; }
	
	REAL x(  ) const
	{ return m_imag.x(  ); }
	REAL y(  ) const
	{ return m_imag.y(  ); }
	REAL z(  ) const
	{ return m_imag.z(  ); }
	REAL w(  ) const
	{ return m_real; }
	
	CQuat & operator*=( const CQuat & q );	
	friend CQuat operator*( const CQuat & q1, const CQuat & q2 );
	friend CPnt operator*( const CQuat & q1, const CPnt & p1 );
	friend CQuat conj( const CQuat & q );
	friend ostream & operator<<( ostream & out, const CQuat & q );
	
	//! CQuat chooseRandom, choose random quaternion
	static CQuat chooseRandom(  );
	
private:
	REAL m_real;		// Real part of the quaterion object
	CPnt m_imag;		// Imaginary part of the quaternion object (3 part)
};  // end CQuat

CPnt randOrient(  );

///////////////////////////////////////////
/////// Inline functions

//!  CPnt -
/*!	Subtracting one xyz coordinate from another  */
inline CPnt 
operator-( const CPnt & c1, const CPnt & c2 )
{
	CPnt c( c1.x( )-c2.x(), c1.y()-c2.y(), c1.z()-c2.z());
	return c;
}

//!  CPnt -
/*! Reversing the direction of a vector */
inline CPnt 
operator-( const CPnt & p1 )
{
	CPnt p( -p1.x( ), -p1.y(), -p1.z()); 
	return p;
}

//!  CPnt +
/*! Adding one xyz coordinate to another */
inline CPnt
operator+( const CPnt & c1, const CPnt & c2 )
{
	CPnt c( c1.x( )+c2.x(), c1.y()+c2.y(), c1.z()+c2.z());
	return c;
}

//!  CPnt * scalar
/*! Multiplying xyz coordinate times a scalar  */
inline CPnt 
operator*( REAL c, const CPnt & p1 )
{
	CPnt p( c * p1.x( ), c * p1.y(), c * p1.z());
	return p;
}

//!  CQuat * CPnt
/*!		Multiplying quaternion by cartesian coordinates, for rotation */
inline CPnt 
operator*( const CQuat & q, const CPnt & r )
{
	CPnt qr = cross( q.m_imag, r );
	return ( r + 2*(q.m_real*qr + cross(q.m_imag, qr )));
}

//!  CQuat *=
/*! Multiplying one quaternion by another */
inline CQuat &
CQuat::operator*=( const CQuat & q )
{
	REAL temp = q.m_real*m_real - dot( q.m_imag, m_imag );
	m_imag =  q.m_real*m_imag + m_real*q.m_imag +	cross( q.m_imag, m_imag );
	m_real = temp;

	return *this;
}

//!  CQuat *
/*!  Multiplying one quaternion by another */
inline CQuat 
operator*( const CQuat & q1, const CQuat & q2 )
{
	CQuat q = q2;
	q *= q1;

	return q;
}

//!  CPnt dot
/*! dot product of two vectors */
inline REAL
dot( const CPnt & p1, const CPnt & p2 )
{
	return ( p1.x( )*p2.x() + p1.y()*p2.y() + p1.z()*p2.z());
}

//!  CPnt cross
/*! Cross product of two vectors  */

inline CPnt 
cross( const CPnt & p1, const CPnt & p2 )
{
	CPnt p( p1.y( )*p2.z() - p1.z()*p2.y(),
	p1.z(  )*p2.x() - p1.x()*p2.z(),
	p1.x(  )*p2.y() - p1.y()*p2.x());

	return p;
}


//!  CPnt torsion
/*! Returning the dihedral angle between 4 points in rad */
inline REAL
torsion( const CPnt & p1, const CPnt & p2, const CPnt & p3, const CPnt & p4 )
{
	return torsion( p1 - p2, p3 - p2, p4 - p3 );
}

//!  CPnt torsion
/*! Returning the dihedral angle between 3 vectors in rad */
inline REAL
torsion( const CPnt & v1, const CPnt & v2, const CPnt & v3 )
{
//      U     W      V
//   o<----o----->o----->o
//

	CPnt A = cross( v1,v2 );
	CPnt B = cross( v3,v2 );

	REAL F = dot( cross(A,B ),v2);
	REAL T = v2.norm(  ) * dot(A,B);

	return atan2( F,T );
}


//!  CPnt <<
/*! Printing out a CPnt object */
inline ostream & 
operator<<( ostream & out, const CPnt & p )
{
	out << "[" << p.x(  ) << "," << p.y() << "," << p.z() <<"]";
	return out;
}

//!  CSpPnt <<
/*! Printing out a CSpPnt object */
inline ostream & 
operator<<( ostream & out, const CSpPnt & p )
{
	out << "[" << p.rho(  ) << "," << p.theta() << "," << p.phi() << "]";
	return out;
}

//!  CSpPnt SphToCart
/*! Converting from spherical to cartesian coordinates */
inline CPnt 
SphToCart( const CSpPnt & s )
{
	REAL x = s.rho(  )*sin(s.theta())*cos(s.phi());
	REAL y = s.rho(  )*sin(s.theta())*sin(s.phi());
	REAL z = s.rho(  )*cos(s.theta());

	return CPnt( x,y,z );
}

//!  CSpPnt CartToSph
/*! Converting from cartesian to spherical coordinates 	*/
inline CSpPnt
CartToSph( const CPnt & c )
{
	REAL theta, phi;
	REAL rho = sqrt( c.x( )*c.x() + c.y()*c.y() + c.z()*c.z());

	if ( rho < fabs(c.z( ))) 
		rho = fabs( c.z( ));

	if ( rho == 0.0 ) 
		theta = 0.0;
	else 
		theta = acos( c.z( )/rho);

	if ( (c.x( ) == 0.0) && (c.y() == 0.0)) 
		phi = 0.0;
	else 
		phi = atan2( c.y( ), c.x());

	return CSpPnt( rho, theta, phi );
}


//!  CQuat conj
/*!		Return complex conjugate of a quaternion */
inline CQuat 
conj( const CQuat & q )
{
	return CQuat( q.m_real, -q.m_imag );
}

//!  CQuat << 
/*! Print out quaternion */
inline ostream & 
operator<<( ostream & out, const CQuat & q )
{
	out << "{" << q.m_real << " " << q.m_imag << "}";
	return out;
}

//!  CQuat chooseRandom 
/*! Choose a random quaternion */
inline CQuat
CQuat::chooseRandom(  )
{
	REAL s = drand48(  );
	REAL sig1 = sqrt( 1.0-s );
	REAL sig2 = sqrt( s );
	REAL t1 = 2.0*drand48(  )*M_PI;
	REAL t2 = 2.0*drand48(  )*M_PI;
	CQuat Q( cos(t2 )*sig2, CPnt(sin(t1)*sig1, cos(t1)*sig1, sin(t2)*sig2));

	return Q;
}

//!  randOrient 
/*!	Choose a random orientation for a CPnt vector */
inline CPnt
randOrient(  )
{
	REAL phi = drand48(  )*2*M_PI;
	REAL u = drand48(  )*2 - 1;

	return CPnt( sqrt(1 - u*u ) * cos(phi), sqrt(1 - u*u) * sin(phi), u);
}

//!  normRand 
/*! Generate a random number  */
inline REAL
normRand(  )
{
	static double saved;
	static bool bSaved = false;
	
	double v1, v2, rsq;
	if ( !bSaved )
	{
		do {
			v1 = 2.0 * drand48(  ) - 1.0;
			v2 = 2.0 * drand48(  ) - 1.0;
			rsq = v1*v1 + v2*v2;
		} while ( rsq >= 1.0 || rsq == 0.0 );
		
		double fac = sqrt( -2.0*log(rsq )/rsq);
		saved = v1*fac;
		bSaved = true;
		return ( REAL ) fac*v2;
	}
	
	else
	{
		bSaved = false;
		return ( REAL ) saved;
	}
}

//!  seedRand 
/*! Generate a random number for random number seed. */
inline void
seedRand( int s )
{
	if ( s == -1 )
	{
		struct timeval tvsd;
		gettimeofday( &tvsd, NULL );
		srand48( tvsd.tv_sec );
	}

	else
		srand48( s );
}

//! resetPos function
/*! reset pos wrt to P
	\param pos a vector of cartesian coordinates to translate
	\param P a cartesian coordinate to translate the positions by */
inline void 
resetPos( vector<CPnt> &pos, CPnt P )
{
	for( int i=0; i < pos.size( ); i++) 
	pos[i] -= P;
}

//! read_timer function
/*! auxiliary function: TIMER */
inline double 
read_timer(   )
{
	static bool initialized = false;
	static struct timeval start, end;
	if(  !initialized  )
	{
		gettimeofday(  &start, NULL  );
		initialized = true;
	}
	gettimeofday(  &end, NULL  );
	return ( end.tv_sec - start.tv_sec ) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

#endif
