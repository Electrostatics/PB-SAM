#include "triexpan.h"

/*###############################################################################
 * #
 * # File: triexpan.cpp
 * #
 * # Date: June 2013
 * #
 * # Description:
 * #
 * # Author: Yap, Felberg
 * #
 * # Copyright ( c )
 * #
 * ################################################################################*/

/*#########################################################*/
/*#########################################################*/
// CTriExpan
/*#########################################################*/
/*#########################################################*/

CTriExpan::CTriExpan( int p, REAL scale )
{
	assert( p >= 0 );
	/*  if ( res > 0 )
	 {
	 m_M[0].reserve( res ); 
	 m_M[1].reserve( res ); 
	 m_M[2].reserve( res ); 
	 }
	 */
	m_M[0].reset( p ); m_M[1].reset(p); m_M[2].reset(p);
	m_M[0].setScale( scale ); m_M[1].setScale(scale);m_M[2].setScale(scale);
}	// end CTriExpan

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

ostream & 
operator<<( ostream & out, const CTriExpan & T )
{
	cout << "\t---000---" << endl;
	cout << T.m_M[0] << endl;
	
	cout << "\t---111---" << endl;
	cout << T.m_M[1] << endl;
	
	cout << "\t---222---" << endl;
	cout << T.m_M[2] << endl;
	
	return out;
}	// end operator<<

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void 
CTriExpan::outputComplex( REAL fact ) const
{
	for( int i=0; i< 3; i++ )
		m_M[i].outputComplex( fact );
}	// end outputComplex

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CTriExpan::rotate( const CQuat & Q, int p )
{
	assert( p <= getRange( ).p2());
	CTriExpan G( p );
	G.setScale(  getScale( ));
	
	for ( int n = 0; n < p; n++ )
	{
		// m = 0 case 
		CPnt ar( m_M[0](n,0 ), m_M[1](n,0), m_M[2](n,0));
		CPnt br = Q*ar;
		G.m_M[0]( n,0 ) = br.x();
		G.m_M[1]( n,0 ) = br.y();
		G.m_M[2]( n,0 ) = br.z();
		
		for ( int m = 1; m <= n; m++ )
		{
			CPnt ar( m_M[0](n,2*m-1 ), m_M[1](n,2*m-1), m_M[2](n,2*m-1));
			CPnt ai( m_M[0](n,2*m ),   m_M[1](n,2*m),   m_M[2](n,2*m));
			CPnt br = Q*ar;
			CPnt bi = Q*ai;
			
			G.m_M[0]( n,2*m-1 ) = br.x();G.m_M[0](n,2*m) =  bi.x();
			G.m_M[1]( n,2*m-1 ) = br.y();G.m_M[1](n,2*m) =  bi.y();
			G.m_M[2]( n,2*m-1 ) = br.z();G.m_M[2](n,2*m) =  bi.z();
		}	// end m
	}	// end n
	this->copy( G, p );
}	 // end rotate

/*#########################################################*/
/*#########################################################*/
// rotate triexpan vectors, assume order 
// has already been incremented elsewhere 
/*#########################################################*/
/*#########################################################*/

void
CTriExpan::incRotate( const CQuat & Q )
{
	int p = getRange(  ).p2();
	assert( p > 0 );
	CTriExpan G( p, getScale( ));
	
	int n = p-1;
	// m = 0 case 
	CPnt ar( m_M[0](n,0 ), m_M[1](n,0), m_M[2](n,0));
	CPnt br = Q*ar;
	G.m_M[0]( n,0 ) = br.x();
	G.m_M[1]( n,0 ) = br.y();
	G.m_M[2]( n,0 ) = br.z();
	
	for ( int m = 1; m <= n; m++ )
	{
		CPnt ar( m_M[0](n,2*m-1 ), m_M[1](n,2*m-1), m_M[2](n,2*m-1));
		CPnt ai( m_M[0](n,2*m ),   m_M[1](n,2*m),   m_M[2](n,2*m));
		CPnt br = Q*ar;
		CPnt bi = Q*ai;
		
		G.m_M[0]( n,2*m-1 ) = br.x(); G.m_M[0](n,2*m) =  bi.x();
		G.m_M[1]( n,2*m-1 ) = br.y(); G.m_M[1](n,2*m) =  bi.y();
		G.m_M[2]( n,2*m-1 ) = br.z(); G.m_M[2](n,2*m) =  bi.z();
	}
	this->copy_p( G, p );
}	// end incRotate

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

CGradExpan::CGradExpan( const CPnt* g, /*const vector<CPnt> &g,*/ 
											 const vector<CPnt> &pos, 
											 int p, bool bKappa, REAL scale, bool bMultipole )
{
	reset( p ); 
	setScale( scale );
	
	if( bMultipole )
		for( int i=0; i<pos.size( ); i++)
		{
			if (  g[i].normsq( ) < 1e-15) continue;
			CSpPnt spos = CartToSph( pos[i] );
			CMExpan M( 1.0, spos, bKappa, p, scale );
			m_M[0] += (  g[i].x( ) * M );
			m_M[1] += (  g[i].y( ) * M );
			m_M[2] += (  g[i].z( ) * M );
		}
	else 
		for( int i=0; i<pos.size( ); i++)
		{
			if (  g[i].normsq( ) < 1e-15) continue;
			CSpPnt spos = CartToSph( pos[i] );
			CLExpan L( 1.0, spos, bKappa, p, scale );
			m_M[0] += (  g[i].x( ) * L );
			m_M[1] += (  g[i].y( ) * L );
			m_M[2] += (  g[i].z( ) * L );
		}
}

