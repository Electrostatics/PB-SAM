#include "rotcoeff.h"

#define EPS_SIN_THETA ( 1e-12 )

REAL CRotCoeff::m_zeta[N_POLES*2][2*N_POLES-1];
REAL CRotCoeff::m_eta[N_POLES*2][4*N_POLES-1];
REAL CRotCoeff::m_Q[2*N_POLES-1][N_POLES][2];

/******************************************************************/
/******************************************************************//**
* Initialize the rotation coefficients (includes coeff equs from 
* 2006 paper
******************************************************************/
void
CRotCoeff::initConstants(  )
{
	for ( int n = 0; n < 2*N_POLES; n++ )
		for ( int m = 0; m <= n; m++ )
			ZETA( n,m ) = sqrt((REAL)(n+m+1)*(n-m+1)/((2*n+1)*(2*n+3)));  // Given as a(n,m) in Lotan, 2006 (eq 1.2)

	for ( int n = 0; n < 2*N_POLES; n++ )
		for ( int m = 0; m <= n; m++ )
		{
			ETA( n,m ) = sqrt((REAL)(n-m-1)*(n-m)/((2*n+1)*(2*n-1)));// Given as b(n,m) in Lotan, 2006 (eq 1.2)

			if ( m != 0 )
				ETA( n,-m ) = -sqrt((REAL)(n+m-1)*(n+m)/((2*n+1)*(2*n-1))); 
		}
	for ( int n = 0; n < 2*N_POLES; n++ )  computeQCoeff();
}	// end initConstants

/******************************************************************/
/******************************************************************//**
*  The rotation coefficient class constructor
******************************************************************/
CRotCoeff::CRotCoeff( bool bGrad ) : m_bSing(false), m_bGrad(bGrad)
{
	m_R.resize( 1 );
	m_R[0].resize( 1 );
	m_R[0][0].resize( 1 );

	if ( m_bGrad )
	{
		m_dR.resize( 1 );
		m_dR[0].resize( 1 );
		m_dR[0][0].resize( 1 );
	}
	m_p = 1;
}	// end CRotCoeff

/******************************************************************/
/******************************************************************//**
*   Allocating space for the R and dR matrices
******************************************************************/
void
CRotCoeff::allocate(  )
{
	int p = m_R.size(  );
	assert(  p+2 <= 2*N_POLES-1 );

	m_R.resize( p+2 );

	m_R[p].resize( 2*p+1 );
	m_R[p+1].resize( 2*p+3 );

	for ( int j = 0; j < 2*p+1; j++ )
		m_R[p][j].resize( p+1 );
	for ( int j = 0; j < 2*p+3; j++ )
		m_R[p+1][j].resize( p+2 );

	if ( m_bGrad )
	{
		m_dR.resize( p+2 );

		m_dR[p].resize( 2*p+1 );
		m_dR[p+1].resize( 2*p+3 );
		for ( int j = 0; j < 2*p+1; j++ )
			m_dR[p][j].resize( p+1 );
		for ( int j = 0; j < 2*p+3; j++ )
			m_dR[p+1][j].resize( p+2 );
	}
}	// end allocate

/******************************************************************/
/******************************************************************//**
*  DeAllocating space for the R and dR matrices
******************************************************************/
void
CRotCoeff::deallocate(  )
{
	int p = m_R.size(  );

	assert( p > 2 );

	for ( int j = 0; j < 2*p-1; j++ )
		m_R[p-1][j].clear(  );
	for ( int j = 0; j < 2*p-3; j++ )
		m_R[p-2][j].clear(  );
	m_R[p-1].clear(  );
	m_R[p-2].clear(  );
	m_R.resize( p-2 );

	if ( m_bGrad )
	{
		for ( int j = 0; j < 2*p-1; j++ )
			m_dR[p-1][j].clear(  );
		for ( int j = 0; j < 2*p-3; j++ )
			m_dR[p-2][j].clear(  );
		m_dR[p-1].clear(  );
		m_dR[p-2].clear(  );
		m_dR.resize( p-2 );
	}
}	// end deallocate

/******************************************************************/
/******************************************************************//**
*  a function to reallocate enough space for the coefficient for 
each pole
******************************************************************/
void
CRotCoeff::reallocate( int p )
{
	assert( p <= N_POLES && p >= 1 );
	if ( p - m_p > 0 )
		for ( int i = m_p; i < p; i++ )
			allocate(  );
	else
		for ( int i = m_p; i > p; i-- )
			deallocate(  );
}	// end reallocate

/******************************************************************/
/******************************************************************//**
*  reset rotation coefficient given an angle for theta, phi and xi
		and a number of poles
		\param theta, a real number of angle theta
		\param phi, a real number for phi angle in rad
		\param xi a
		\param p an integer indicating the number of poles 
******************************************************************/
void 
CRotCoeff::reset( REAL theta, REAL phi, REAL xi, int p )
{
	assert( p > 0 && p <= N_POLES );
	initParams( theta, phi, xi );
	reallocate( p );

	m_SH = CSHExpan( m_theta, m_phi, 2*p-1, true );

	m_p = p;  

	if ( !isSingular( ))
	{
		computeCoeff(  );
		if ( m_bGrad ) computeGradCoeff();
	}
}	// end reset

/******************************************************************/
/******************************************************************//**
*  reset: Resets the rotation coefficients for a given number of poles
	\param Q a quaternion for rotation
	\param p an int describing the number of poles
******************************************************************/
void
CRotCoeff::reset( const CQuat & Q, int p )
{
	assert( p > 0 && p <= N_POLES );

	REAL theta = acos( 1.0 - 2*Q.x( )*Q.x() - 2*Q.y()*Q.y());

	REAL phi, xi;
	if ( fabs(theta ) > 1e-12) 
	{ 
		phi = atan2( Q.y( )*Q.z() + Q.w()*Q.x(), Q.x()*Q.z() - Q.w()*Q.y());
		xi = atan2( Q.y( )*Q.z() - Q.w()*Q.x(), Q.x()*Q.z() + Q.w()*Q.y());
	}
	else
	{
		theta = 0.0;
		phi = 0.0; //phi = M_PI/2.0;
		xi = 2*acos( Q.w( ));   //xi = -M_PI/2.0;
	}
	reset( theta, phi, xi, p );
}	// end reset

/******************************************************************/
/******************************************************************//**
*  Initializing parameters in RotCoeff object
******************************************************************/
void
CRotCoeff::initParams( REAL theta, REAL phi, REAL xi )
{
	//cout <<"rotcoeff initparams theta "<<theta
	//<<" phi "<<phi<<"  xi "<<xi<<endl;
	m_theta = theta; m_phi = phi; m_xi = xi;

	m_sint = sin( m_theta );
	if ( m_sint < EPS_SIN_THETA )
	{
		m_bSing = true;
		if ( m_theta < M_PI/2 )
		{
			m_theta = 0.0;
			m_phi = 0.0;
			m_sint = 0.0;
			m_cost = 1.0;
			m_exphi = Complex( 1.0,0.0 );
		}
		else
		{
			m_theta = M_PI;
			m_phi = 0.0;
			m_sint = 0.0;
			m_cost = -1.0;
			m_exphi = Complex( 1.0,0.0 );
		}
	}
	else
	{
		m_bSing = false;
		m_cost = cos( m_theta );
		m_cott = m_cost/m_sint;
		m_exphi = Complex( cos(phi ),sin(phi));

		//check for near zero values - enghui
		/*
		if( m_cost < EPS_SIN_THETA  ) {
		m_cost = 0.0;
		m_sint = 1.0;
		m_cott = 0.0;
		}
		if( cos(phi ) < EPS_SIN_THETA) m_exphi = Complex(0.0,sin(phi));
		if( sin(phi ) < EPS_SIN_THETA) m_exphi = Complex(cos(phi),0.0);
		*/
	}
	m_exiphi = conj( m_exphi );
	m_exxi = Complex( cos(xi ),sin(xi));	  
}	// end initParams

/*#########################################################*/
/*#########################################################*/
// !!! Note: used only for general rotation 
// e.g. molecular rotation. 
// Uses xi from reset(  ) for singular cases.
// Apply the rotation operator to the MP coeffs
/*#########################################################*/
/*#########################################################*/

void
CRotCoeff::rotateWithXi( const CExpan & Min, CExpan & Mout, int p1, int p2, 
			bool bFor ) const 
{
	assert( Min.getRange( ).p2() >= p2 && p2 <= N_POLES && p2 <= m_p);
	assert( p1 == 1 || Mout.getRange( ).p2() == p1 - 1);
	assert( p1 != 0 && p1 <= p2 );

	Mout.setScale( Min.getScale( ));
	Mout.setRange( p2 );

	if ( isSingular( ))
	{
		REAL xi = ( bFor ? m_xi : -m_xi );

		if ( m_cost == -1.0 )
		{
			for ( int n = p1-1; n < p2; n++ )
			{
				Mout( n,0 ) = (n % 2 == 0 ? Min(n,0) : -Min(n,0));

				int s = ( n % 2 == 0 ? 1 : -1 );
				for ( int m = 1; m <= n; m++ ) 
				{
					Complex incplx(  Min(n,2*m-1 ) , Min(n,2*m) );
					Complex ex_imxi(  cos(m*xi ), sin(m*xi)); 
					Complex outcplx = ex_imxi * incplx;
					Mout( n,2*m-1 ) = s*outcplx.real();
					Mout( n,2*m )   = -s*outcplx.imag();
				}	// end m
			}	// end n
		}	// end if m_cost
		else
		{
			for ( int n = p1-1; n < p2; n++ )
			{
				Mout( n,0 ) = Min(n,0);		  

				for ( int m = 1; m <= n; m++ ) 
				{
					Complex incplx(  Min(n,2*m-1 ) , Min(n,2*m) );
					Complex ex_imxi(  cos(m*xi ), sin(m*xi)); 
					Complex outcplx = ex_imxi * incplx;
					Mout( n,2*m-1 ) = outcplx.real();
					Mout( n,2*m )   = outcplx.imag();
				} // end m
			}// end n
		} // end else
	} // end isSingular
	else // Not Singular
	{
		if ( bFor )
		{
			for ( int n = p1 - 1; n < p2; n++ )
				for ( int m = 0; m <= n; m++ )
				{
					Complex incplx = Complex( Min(n,0 ), 0.0);
					Complex outcplx = CExpan::SQRT2 * ROT( n,0,m ) * incplx;

					for ( int s = 1; s <= n; s++ ) 
					{
						incplx = Complex(  Min(n,2*s-1 ) , Min(n,2*s) );
						outcplx += ( ROT(n,s,m )*incplx + ROT(n,-s,m)*conj(incplx)); 
					}
					if ( m == 0 )
						Mout( n,0 ) = outcplx.real() * CExpan::ISQRT2;		  
					else 
					{
						Mout( n,2*m-1 ) = outcplx.real();
						Mout( n,2*m )   = outcplx.imag();
					}
				}	// end n m
		}
		else // TRANSPOSED
		{
			for ( int n = p1 - 1; n < p2; n++ )
				for ( int s = 0; s <= n ; s++ )
				{
					Complex incplx = Complex( Min(n,0 ), 0.0);
					Complex outcplx = CExpan::SQRT2 * ROT( n,-s,0 ) * incplx;

					for ( int m = 1; m <= n; m++ ) 
					{
						incplx = Complex(  Min(n,2*m-1 ) , Min(n,2*m) );
						outcplx += conj( ROT(n,s,m ))*incplx + ROT(n,-s,m)*conj(incplx);
					}

					if ( s == 0 ) Mout(n,0) = outcplx.real() * CExpan::ISQRT2;
					else 
					{
						Mout( n,2*s-1 ) = outcplx.real();
						Mout( n,2*s )   = outcplx.imag();
					}
				} // end n and s
	}	// end else (  transposed  )
}	// end else not singular
} // end rotateWithXi

/*#########################################################*/
/*#########################################################*/
// !!! Note: used only for rotation as 
// part of RSR transform, 
// Has inbuild assumptions about values of 
// xi for singular cases.
// Apply the rotation operator to the MP coeffs
/*#########################################################*/
/*#########################################################*/

void
CRotCoeff::rotate( const CExpan & Min, CExpan & Mout, int p1, int p2, 
			bool bFor ) const
{
	assert( Min.getRange( ).p2() >= p2 && p2 <= N_POLES && p2 <= m_p);
	assert( p1 == 1 || Mout.getRange( ).p2() == p1 - 1);
	assert( p1 != 0 && p1 <= p2 );

	Mout.setScale( Min.getScale( ));
	Mout.setRange( p2 );

	if ( isSingular( ))
	{
		if ( m_cost == -1.0 )
		{
			//	  cout<<"in singular case pi"<<endl;
			for ( int n = p1 - 1; n < p2; n++ ) 
			{
				Mout( n,0 ) = (n % 2 == 0 ? Min(n,0) : -Min(n,0));

				int s = ( n % 2 == 0 ? 1 : -1 );
				for ( int m = 1; m < 2*n+1; m++, s = -s ) 
				Mout( n,m ) = s*Min(n,m);
			}	// end n
		}	// end if m_cost
		else
		{
			for ( int n = p1-1; n < p2; n++ )
				for ( int m = 0; m < 2*n+1; m++ ) 
				{
					Mout( n,m ) = Min(n,m);
				}
		}	// end else
	}	// end isSingular
	else // Not Singular
	{
		if ( bFor )
		{
			for ( int n = p1 - 1; n < p2; n++ )
				for ( int m = 0; m <= n; m++ )
				{
					Complex incplx = Complex( Min(n,0 ), 0.0);
					Complex outcplx = ROT( n,0,m ) * incplx;
					for ( int s = 1; s <= n; s++ ) 
					{
						incplx = Complex(  Min(n,2*s-1 ) , Min(n,2*s) );
						outcplx += ( ROT(n,s,m )*incplx + ROT(n,-s,m)*conj(incplx)); 
					}
					if ( m == 0 )
						Mout( n,0 ) = outcplx.real();		  
					else 
					{
						Mout( n,2*m-1 ) = outcplx.real();
						Mout( n,2*m )   = outcplx.imag();
					}
				}	// end n m
		}	// end if Bfor
		else // TRANSPOSED
		{
			for ( int n = p1 - 1; n < p2; n++ )
				for ( int s = 0; s <= n ; s++ )
				{
					Complex incplx( Min(n,0 ), 0.0);
					Complex outcplx = CExpan::SQRT2 * ROT( n,-s,0 ) * incplx;

					for ( int m = 1; m <= n; m++ ) 
					{
						incplx = Complex(  Min(n,2*m-1 ) , Min(n,2*m) );
						outcplx += conj( ROT(n,s,m ))*incplx + ROT(n,-s,m)*conj(incplx);
					}
					if ( s == 0 ) Mout(n,0) = outcplx.real() * CExpan::ISQRT2;
					else 
					{
						Mout( n,2*s-1 ) = outcplx.real();
						Mout( n,2*s )   = outcplx.imag();
					}
				}	// end for ns
		}	// end else transposed
	}	// end else not singular
}	// end rotate

/******************************************************************/
/******************************************************************//**
* Apply the derivative of the rotation operator with respect to THETA
    to the MP coeffs
******************************************************************/
void
CRotCoeff::dRotateT( const CExpan & Min, CExpan & Mout, int p1, int p2,
	bool bFor ) const
{
	if ( !m_bGrad )
	{
		cout << "dRdt: This operator has no derivatives!!!" << endl;
		return;
	}

	assert( Min.getRange( ).p2() >= p2 && p2 <= N_POLES && p2 <= m_p);
	assert( p1 == 1 || Mout.getRange( ).p2() == p1 - 1);
	assert( p1 != 0 && p1 <= p2 );

	Mout.setRange( p2 );
	Mout.setScale( Min.getScale( ));

	if ( isSingular( ))
	{
		dRotateTSing( Min, Mout, p1, p2 );
		if ( !bFor )
			Mout.recip(  );
	}
	else
	{
		if ( bFor )
		{
			for ( int n = p1 - 1; n < p2; n++ )
				for ( int m = 0; m <= n; m++ )
				{
					Complex incplx(  Min(n,0 ), 0.0);
					Complex outcplx = CExpan::SQRT2 * dROT( n,0,m ) * incplx;

					for ( int s = 1; s <= n; s++ )
					{ 
						incplx = Complex(  Min(n,2*s-1 ) , Min(n,2*s) );
						outcplx += ( dROT(n,s,m )*incplx + dROT(n,-s,m)*conj(incplx)); 
						//Mout( n,m ) += (dROT(n,s,m)*Min(n,s) + dROT(n,-s,m)*Min(n,-s));
					}
					if ( m == 0 )
						Mout( n,0 ) = outcplx.real() * CExpan::ISQRT2;		  
					else 
					{
						Mout( n,2*m-1 ) = outcplx.real();
						Mout( n,2*m )   = outcplx.imag();
					}
				}	// end m n
		}	// end bFor
		else // TRANSPOSED
		{
			for ( int n =  p1 - 1; n < p2; n++ )
				for ( int s = 0; s <= n; s++ )
				{
					Complex incplx(  Min(n,0 ), 0.0);
					Complex outcplx = CExpan::SQRT2 * dROT( n,-s,0 ) * incplx;
					//Mout( n,s ) = dROT(n,-s,0)*Min(n,0);

					for ( int m = 1; m <= n; m++ )
					{
						incplx = Complex(  Min(n,2*m-1 ) , Min(n,2*m) );
						outcplx += conj( dROT(n,s,m ))*incplx + dROT(n,-s,m)*conj(incplx);
						/*Mout( n,s ) += conj(dROT(n,s,m))*Min(n,m) +
									dROT( n,-s,m )*Min(n,-m);*/
					}
					if ( s == 0 ) Mout(n,0) = outcplx.real() * CExpan::ISQRT2;
					else 
					{
						Mout( n,2*s-1 ) = outcplx.real();
						Mout( n,2*s )   = outcplx.imag();
					}		  
				}	// end n and s
		}	// end else transposed
	}	// end else not singular
}	// end dRotateT

/******************************************************************/
/******************************************************************//**
*  Apply the derivative of the rotation operator with respect to PHI 
			to the MP coeffs
******************************************************************/
void
CRotCoeff::dRotateP( const CExpan & Min, CExpan & Mout, int p1, int p2,
			bool bFor ) const
{
	if ( !m_bGrad )
	{
		cout << "dRdp: This operator has no derivatives!!!" << endl;
		return;
	}

	assert( Min.getRange( ).p2() >= p2 && p2 <= N_POLES && p2 <= m_p);
	assert( p1 == 1 || Mout.getRange( ).p2() == p1 - 1);
	assert( p1 != 0 && p1 <= p2 );

	Mout.setRange( p2 );
	Mout.setScale( Min.getScale( ));

	if ( isSingular( ))
	{
		dRotatePSing( Min, Mout, p1, p2 );
		if ( !bFor && m_cost == 1.0 )
			Mout.recip(  );
	}
	else
	{
		if ( bFor )
		{
			for ( int n = p1 - 1; n < p2; n++ )
				for ( int m = 0; m <= n; m++ )
				{
					Complex incplx;
					Complex outcplx( 0.0,0.0 ); // s = 0

					for ( int s = 1; s <= n; s++ )
					{
						incplx = Complex(  Min(n,2*s-1 ) , Min(n,2*s) );
						outcplx -=  ( ROT(n,s,m )*derv(incplx,s) + 
						ROT( n,-s,m )*derv(conj(incplx),-s)); 
					}

					if ( m == 0 )
						Mout( n,0 ) = outcplx.real() * CExpan::ISQRT2;		  
					else 
					{
						Mout( n,2*m-1 ) = outcplx.real();
						Mout( n,2*m )   = outcplx.imag();
					}
				}
		}
		else
		{
			for ( int n = p1 - 1 ; n < p2; n++ )
				for ( int s = 0; s <= n; s++ )
				{
					Complex incplx = Complex( Min(n,0 ),0);
					Complex outcplx = CExpan::SQRT2 * ROT( n,-s,0 )
					*derv( incplx,s ) ; // m = 0

					for ( int m = 1; m <= n; m++ )
					{
						incplx = Complex(  Min(n,2*m-1 ) , Min(n,2*m) );
						outcplx += conj( ROT(n,s,m ))*derv(incplx,s) + 
						ROT( n,-s,m )*derv(conj(incplx),s);
					}

				if ( s == 0 ) Mout(n,0) = outcplx.real() * CExpan::ISQRT2;
				else 
				{
					Mout( n,2*s-1 ) = outcplx.real();
					Mout( n,2*s )   = outcplx.imag();
				}	
			}	// end s n
		}	// end else
	}  // end  else
}	// end dRotateP

/******************************************************************/
/******************************************************************//**
*  Apply the derivative of the rotation operator with respect to PHI 
			to the MP coeffs in the singular case sin(theta) = 0.0, 
			given as equation 1.11 in paper Lotan 2006.
******************************************************************/
void 
CRotCoeff::dRotateTSing( const CExpan & Min, 
			CExpan & Mout, int p1, int p2 ) const
{
	assert( Min.getRange( ).p2() >= p2 && p2 <= N_POLES && p2 <= m_p);
	assert( p1 == 1 || Mout.getRange( ).p2() == p1 - 1);  
	assert( p1 != 0 && p1 <= p2 );

	if ( p1 == 1 )
	{
		Mout( 0,0 ) = 0.0;
		p1++;
	}

	if ( m_cost == 1.0 )
	{
		for ( int n = p1 - 1; n < p2; n++ )
		{
			Mout( n,0 ) = 2.0*m_Q[n][0][1]*Min(n,1) * CExpan::ISQRT2;
			for ( int m = 1; m < n; m++ )
			{
				Complex incplx0, outcplx;
				if( m == 1 ) incplx0 = Complex(Min(n,0), 0.0) * CExpan::SQRT2; 
				else incplx0 = Complex( Min(n,2*(m-1 )-1), Min(n,2*(m-1)));
				Complex incplx1( Min(n,2*(m+1 )-1), Min(n,2*(m+1)));

				outcplx = m_Q[n][m][0]*incplx0 + m_Q[n][m][1]*incplx1; 
				Mout( n,2*m-1 ) = outcplx.real();
				Mout( n,2*m   ) = outcplx.imag();
			}

			Complex incplx0;
			if( n == 1 ) incplx0 = Complex(Min(n,0), 0.0) * CExpan::SQRT2; 
			else incplx0 = Complex( Min(n,2*(n-1 )-1), Min(n,2*(n-1)));
			Complex outcplx = m_Q[n][n][0]*incplx0; 
			Mout( n,2*n-1 )= outcplx.real();
			Mout( n,2*n ) = outcplx.imag();
		} // end n
	}	// end if m_cost
	else
	{
		REAL s = -1.0;
		for ( int n = p1 - 1; n < p2; n++, s = -s )
		{
			Mout( n,0 ) = 2.0*s*m_Q[n][0][1]*Min(n,1) * CExpan::ISQRT2;

			for ( int m = 1; m < n; m++ )
			{
				Complex incplx0,outcplx;
				if( m == 1 ) incplx0 = Complex(Min(n,0), 0.0) * CExpan::SQRT2; 
				else incplx0 = Complex( Min(n,2*(m-1 )-1), Min(n,2*(m-1)));
				Complex incplx1( Min(n,2*(m+1 )-1), Min(n,2*(m+1)));

				outcplx =  s*( m_Q[n][m][0] * conj(incplx0 ) +
				m_Q[n][m][1] * conj( incplx1 ) );
				Mout( n,2*m-1 ) = outcplx.real();
				Mout( n,2*m ) = outcplx.imag();
			}	// end for m
			Complex incplx0;
			if( n == 1 ) incplx0 = Complex(Min(n,0), 0.0) * CExpan::SQRT2;
			else incplx0 = Complex( Min(n,2*(n-1 )-1), Min(n,2*(n-1)));
			Complex outcplx = s*m_Q[n][n][0]* conj( incplx0 );
			Mout( n,2*n-1 ) = outcplx.real();
			Mout( n,2*n ) = outcplx.imag();
		}	// end for n
	}	// end else
}	// end dRotateTSing


/******************************************************************/
/******************************************************************//**
*  Apply the derivative of the rotation operator with respect to PHI 
			to the MP coeffs in the singular case sin(theta) = 0.0
 !!! Only used as part of RSR transform. 
 Singular cases has inbuild assumption about values of xi
 Apply the derivative of the rotation operator with respect to PHI 
 to the MP coeffs in the singular case sin( theta ) = 0.0
******************************************************************/
void 
CRotCoeff::dRotatePSing( const CExpan & Min, 
			CExpan & Mout, int p1, int p2 ) const
{
	assert( Min.getRange( ).p2() >= p2 && p2 <= N_POLES && p2 <= m_p);
	assert( p1 == 1 || Mout.getRange( ).p2() == p1 - 1); 
	assert( p1 != 0 && p1 <= p2 );

	if ( p1 == 1 )
	{
		Mout( 0,0 ) = 0.0;
		p1++;
	}

	if ( m_cost == 1.0 )
	{
		for ( int n = p1 - 1; n < p2; n++ )
		{
			Mout( n,0 ) = 2.0*m_Q[n][0][1]*Min(n,2) * CExpan::ISQRT2; 
			//	  Mout( n,0 ) = Complex(2.0*m_Q[n][0][1]*Min(n,1).imag(), 0.0); 
			for ( int m = 1; m < n; m++ )
			{
				Complex incplx0, outcplx;
				if( m == 1 ) incplx0 = Complex(Min(n,0), 0.0) * CExpan::SQRT2; 
				else incplx0 = Complex( Min(n,2*(m-1 )-1), Min(n,2*(m-1)));
				Complex incplx1( Min(n,2*(m+1 )-1), Min(n,2*(m+1)));

				outcplx = Complex( 0.0,m_Q[n][m][0] )* incplx0 - 
				Complex( 0.0,m_Q[n][m][1] )*incplx1; 

				Mout( n,2*m-1 ) = outcplx.real();
				Mout( n,2*m ) = outcplx.imag();
			}

			Complex incplx0;
			if( n == 1 ) incplx0 = Complex(Min(n,0), 0.0) * CExpan::SQRT2; 
			else incplx0 = Complex( Min(n,2*(n-1 )-1), Min(n,2*(n-1)));

			Complex outcplx = Complex( 0.0, m_Q[n][n][0] )*incplx0; 
			Mout( n,2*n-1 )= outcplx.real();
			Mout( n,2*n ) = outcplx.imag();
		} // end n
	}
	else
	{
		REAL s = 1.0;
		for ( int n = p1 - 1; n < p2; n++, s = -s )
		{
			Mout( n,0 ) = s*2*m_Q[n][0][1]*Min(n,2) * CExpan::ISQRT2;

			for ( int m = 1; m < n; m++ )
			{
				Complex incplx0,outcplx;
				if( m == 1 ) incplx0 = Complex(Min(n,0), 0.0) * CExpan::SQRT2; 
				else incplx0 = Complex( Min(n,2*(m-1 )-1), Min(n,2*(m-1)));
				Complex incplx1( Min(n,2*(m+1 )-1), Min(n,2*(m+1)));

				outcplx =  s*( Complex(0.0,m_Q[n][m][1] )* conj(incplx1) -
				Complex( 0.0,m_Q[n][m][0] )* conj(incplx0));
				Mout( n,2*m-1 ) = outcplx.real();
				Mout( n,2*m ) = outcplx.imag();
			}

			Complex incplx0;
			if( n == 1 ) incplx0 = Complex(Min(n,0), 0.0) * CExpan::SQRT2; 
			else incplx0 = Complex( Min(n,2*(n-1 )-1), Min(n,2*(n-1)));
			Complex outcplx = Complex( 0.0,-s*m_Q[n][n][0] )* conj(incplx0);
			Mout( n,2*n-1 ) = outcplx.real();
			Mout( n,2*n ) = outcplx.imag();
		} // end n
	} // end else
}	// end dRotatePSing

/******************************************************************/
/******************************************************************//**
*   Function to compute rotation coefficients. Using components from
appendix of Lotan 2006 paper. 
******************************************************************/
void
CRotCoeff::computeCoeff(  )
{
	Complex dum1, dum2, dum3; 

	ROT( 0,0,0 ) = m_SH.comp(0,0); // access SH as a complex number

	for ( int n = 1; n < 2*m_p-1; n++ )
	{
		ROT( n,0,0 ) = m_SH.comp(n,0);
		for ( int s = 1; s <= n; s++ )
		{	 
			ROT( n,s,0 ) = m_SH.comp(n,-s);
			ROT( n,-s,0 ) = m_SH.comp(n,s);	  
		}
	}

	m_r1 = -0.5*( 1 + m_cost );
	m_r2 = 0.5*( 1 - m_cost );
	m_r3 = -m_sint;

	for ( int m = 0; m < m_p-1; m++ )
		for ( int n = m+2; n < 2*m_p-m-1; n++ )
			for ( int s = -n+1; s < n; s++ )
			{	
				dum1 = ( m_r1*ETA(n,s-1 ))*ROT(n,s-1,m);
				dum2 = ( m_r2*ETA(n,-s-1 ))*ROT(n,s+1,m);
				dum3 = ( m_r3*ZETA(n-1,abs(s )))*ROT(n,s,m);
				ROT( n-1,s,m+1 ) = m_exxi*(dum1*m_exiphi+dum2
									*m_exphi+dum3 )/ETA( n,m);
			}
}	// end computeCoeff

/******************************************************************/
/******************************************************************//**
*     Function to compute gradient of rotation coefficients 
******************************************************************/
void
CRotCoeff::computeGradCoeff(  )
{
	Complex dum1, dum2, dum3;

	dROT( 0,0,0 ) = Complex(0.0,0.0);
	for ( int n = 1; n < 2*m_p-1; n++ )
	{
		dROT( n,0,0 ) = -sqrt((REAL)n*(n+1))*m_exiphi*m_SH.comp(n,1);
		for ( int s = 1; s <  n; s++ )
		{	 
			dROT( n,-s,0 ) = (s*m_cott)*m_SH.comp(n,s) -
			sqrt( (REAL )(n-s)*(n+s+1))*m_exiphi*m_SH.comp(n,s+1);
			dROT( n,s,0 ) = conj(dROT(n,-s,0));
		}

		dROT( n,-n,0 ) = (n*m_cott)*m_SH.comp(n,n);
		dROT( n,n,0 ) = conj(dROT(n,-n,0));
	}

	m_dr1 = 0.5*m_sint;
	m_dr2 = m_dr1;
	m_dr3 = -m_cost;
	for ( int m = 0; m < m_p-1; m++ )
		for ( int n = m+2; n < 2*m_p-m-1; n++ )
			for ( int s = -n+1; s < n; s++ )
			{
				dum1 = ETA( n,s-1 )*(m_r1*dROT(n,s-1,m) + m_dr1*ROT(n,s-1,m));
				dum2 = ETA( n,-s-1 )*(m_r2*dROT(n,s+1,m) + m_dr2*ROT(n,s+1,m));
				dum3 = ZETA( n-1,abs(s ))*(m_r3*dROT(n,s,m) + m_dr3*ROT(n,s,m));
				dROT( n-1,s,m+1 ) = m_exxi*(dum1*m_exiphi
							+dum2*m_exphi+dum3 )/ETA( n,m);
			}
}	// end computeGradCoeff

/******************************************************************/
/******************************************************************//**
*  Function to compute additional rotational coefficients when the 
			number of poles is increased. 
******************************************************************/
void
CRotCoeff::computeIncCoeff(  )
{
	Complex dum1, dum2, dum3;

	for ( int n = 2*m_p-3; n < 2*m_p-1; n++ )
	{
		ROT( n,0,0 ) = m_SH.comp(n,0);
		for ( int s = 1; s < n; s++ )
		{	 
			ROT( n,s,0 ) = m_SH.comp(n,-s);
			ROT( n,-s,0 ) = m_SH.comp(n,s);
		}
		ROT( n,n,0 ) = m_SH.comp(n,-n);
		ROT( n,-n,0 ) = m_SH.comp(n,n);
	}

	for ( int m = 0; m < m_p-2; m++ )
		for ( int n = 2*m_p-m-3; n < 2*m_p-m-1; n++ )
			for ( int s = -n+1; s < n; s++ )
			{
				dum1 = ( m_r1*ETA(n,s-1 ))*ROT(n,s-1,m);
				dum2 = ( m_r2*ETA(n,-s-1 ))*ROT(n,s+1,m);
				dum3 = ( m_r3*ZETA(n-1,abs(s )))*ROT(n,s,m);
				ROT( n-1,s,m+1 ) = m_exxi*(dum1*m_exiphi
								+dum2*m_exphi+dum3 )/ETA( n,m);
			}

	int m = m_p-2;
	for ( int n = m+2; n < 2*m_p-m-1; n++ )
		for ( int s = -n+1; s < n; s++ )
		{
			dum1 = ( m_r1*ETA(n,s-1 ))*ROT(n,s-1,m);
			dum2 = ( m_r2*ETA(n,-s-1 ))*ROT(n,s+1,m);
			dum3 = ( m_r3*ZETA(n-1,abs(s )))*ROT(n,s,m);
			ROT( n-1,s,m+1 ) = m_exxi*(dum1*m_exiphi 
								+ dum2*m_exphi + dum3 )/ETA( n,m);
		}
}	// end computeIncCoeff

/******************************************************************/
/******************************************************************//**
* Function to compute additional gradient rot coefficients when the 
			number of poles is increased.
******************************************************************/
void
CRotCoeff::computeIncGradCoeff(  )
{ 
	Complex dum1, dum2, dum3;

	for ( int n = 2*m_p-3; n < 2*m_p-1; n++ )
	{
		dROT( n,0,0 ) = -sqrt((REAL)n*(n+1))*m_exiphi*m_SH.comp(n,1);
		for ( int s = 1; s < n; s++ )
		{	 
			dROT( n,-s,0 ) = (s*m_cott)*m_SH.comp(n,s) -
			sqrt( (REAL )(n-s)*(n+s+1))*m_exiphi*m_SH.comp(n,s+1);
			dROT( n,s,0 ) = conj(dROT(n,-s,0));
		}

		dROT( n,-n,0 ) = (n*m_cott)*m_SH.comp(n,n);
		dROT( n,n,0 ) = conj(dROT(n,-n,0));
	}

	for ( int m = 0; m < m_p-2; m++ )
		for ( int n = 2*m_p-m-3; n < 2*m_p-m-1; n++ )
			for ( int s = -n+1; s < n; s++ )
			{
				dum1 = ETA( n,s-1 )*(m_r1*dROT(n,s-1,m) + m_dr1*ROT(n,s-1,m));
				dum2 = ETA( n,-s-1 )*(m_r2*dROT(n,s+1,m) + m_dr2*ROT(n,s+1,m));
				dum3 = ZETA( n-1,abs(s ))*(m_r3*dROT(n,s,m) + m_dr3*ROT(n,s,m));
				dROT( n-1,s,m+1 ) = m_exxi*(dum1*m_exiphi
							+dum2*m_exphi+dum3 )/ETA( n,m);
			}

	int m = m_p-2;
	for ( int n = m+2; n < 2*m_p-m-1; n++ )
		for ( int s = -n+1; s < n; s++ )
		{
			dum1 = ETA( n,s-1 )*(m_r1*dROT(n,s-1,m) + m_dr1*ROT(n,s-1,m));
			dum2 = ETA( n,-s-1 )*(m_r2*dROT(n,s+1,m) + m_dr2*ROT(n,s+1,m));
			dum3 = ZETA( n-1,abs(s ))*(m_r3*dROT(n,s,m) + m_dr3*ROT(n,s,m));
			dROT( n-1,s,m+1 ) = m_exxi*(dum1*m_exiphi
							+dum2*m_exphi+dum3 )/ETA( n,m);
		}
}	// end computeIncGradCoeff

/******************************************************************/
/******************************************************************//**
*  Computing Q coefficients, which are only used for computing 
the dR when sint -> 0
******************************************************************/
void
CRotCoeff::computeQCoeff(  )
{
	vector<REAL> SH( 2*N_POLES-1 );
	CSHExpan::specialSH( SH, 2*N_POLES-1, 1.0 );

	REAL temp;
	for ( int n = 1; n < 2*N_POLES-1; n++ )
	{
		m_Q[n][0][0] = 0.0;
		m_Q[n][0][1] = SH[n];
	}

	for ( int n = 2; n < 2*N_POLES-1; n++ )
	{
		m_Q[n-1][1][0] = ( ETA(n,-1 )*m_Q[n][0][1] 
							+  ZETA( n-1,0 ))/ETA(n,0);
		if ( n > 2 )
			m_Q[n-1][1][1] = ( ETA(n,1 )/ETA(n,0))*m_Q[n][0][1];
		else
			m_Q[n-1][1][1] = 0.0;
	}

	for ( int m = 1; m < N_POLES-1; m++ )
	for ( int n = m+2; n < 2*N_POLES-m-1; n++ )
	{
		m_Q[n-1][m+1][0] = ( ETA(n,m-1 )*m_Q[n][m][0] 
							+  ZETA( n-1,m ))/ETA(n,m);
		if ( n > m+2 )
			m_Q[n-1][m+1][1] =  ( ETA(n,m+1 )/ETA(n,m))*m_Q[n][m][1];
		else
			m_Q[n-1][m+1][1] = 0.0;
	}
}	// end computeQCoeff

/******************************************************************/
/******************************************************************//**
*  Function to increase the number of poles and recompute rotation
coefficients
******************************************************************/
void
CRotCoeff::incOrder(  )
{
	allocate(  );
	m_p++;

	if ( !isSingular( ))
	{
		while ( m_SH.getRange( ).p2() < 2*m_p-1)
		m_SH.inc(  );

		computeIncCoeff(  );
		if ( m_bGrad )	computeIncGradCoeff();
	}
}	// end incOrder

/******************************************************************/
/******************************************************************//**
*    Print out rotation coefficients
******************************************************************/
void
CRotCoeff::outputRot( int p ) const
{
	if ( p > m_p )
		p = m_p;

	cout << "****ROT COEFFICIENTS****" << endl;
	for ( int m = 0; m < m_p; m++ )
	{
		cout << "\t---m = " << m << "---" << endl;
		for ( int n = m; n < m_p; n++ )
		{
			for ( int s = -n; s <= n; s++ )
			{
				REAL r = fabs( ROT(n,s,m ).real())>1e-15 ? 
				ROT( n,s,m ).real() : 0;
				REAL im = fabs( ROT(n,s,m ).imag())>1e-15 ? 
				ROT( n,s,m ).imag() : 0;
				cout << "( " << r << "," << im << " ) | ";
			}
			cout << endl;
		}
	}
}	// end outputRot

/******************************************************************/
/******************************************************************//**
*    Print out derivatives of rot coefficients
******************************************************************/
void
CRotCoeff::outputdRot( int p ) const
{
	if ( p > m_p )
		p = m_p;

	cout << "****dROT COEFFICIENTS****" << endl;
	for ( int m = 0; m < m_p; m++ )
	{
		cout << "\t---m = " << m << "---" << endl;
		for ( int n = m; n < m_p; n++ )
		{
			for ( int s = -n; s <= n; s++ )
			{
				REAL r = fabs( dROT(n,s,m ).real())>1e-15 ? 
				dROT( n,s,m ).real() : 0;
				REAL im = fabs( dROT(n,s,m ).imag())>1e-15 ? 
				dROT( n,s,m ).imag() : 0;
				cout << "( " << r << "," << im << " ) | ";
			}
			cout << endl;
		}
	}
}	// end outputdRot

/******************************************************************//**
* Function to check the rotation matrix size and ensure that there is
a matrix coefficient at a given position.
s_ = actual index address ( should be = n+s when called )
******************************************************************/
bool
CRotCoeff::checkR( int n, int s, int m ) const
{
	int sN = m_R.size(  );
	bool bsN = ( n < sN  );
	if( !bsN )
	{
		cout <<"R too small: n="<<n<<" "<<sN<<endl; 
		exit( 0 ); 
	}

	int sS = m_R[n].size(  );
	bool bsS=(  s < sS );
	if( !bsS )
	{
		cout <<"R too small: s= "<<s<<" "<<sS<<endl; 
		exit( 0 ); 
	}

	int sM = m_R[n][s].size(  );
	bool bsM= ( m < sM  );
	if( !bsM )
	{
		cout <<"R too small: m="<<m<<" "<<sM<<endl; 
		exit( 0 ); 
	}
	return true;
}// end checkR

/******************************************************************//**
* Function to check the derivative rotation matrix size and 
ensure that there is a matrix coefficient at a given position.
s_ = actual index address ( should be = n+s when called )
******************************************************************/
bool
CRotCoeff::checkdR( int n, int s, int m ) const
{
	int sN = m_dR.size(  );
	bool bsN = ( n < sN  );
	if( !bsN )
	{
		cout <<"dR too small: n="<<n<<" "<<sN<<endl; 
		exit( 0 ); 
	}

	int sS = m_dR[n].size(  );
	bool bsS=(  s < sS  );
	if( !bsS )
	{
		cout <<"dR too small: s= "<<s<<" "<<sS<<endl; 
		exit( 0 );
	}

	int sM = m_dR[n][s].size(  );
	bool bsM= ( m < sM  );
	if( !bsM )
	{
		cout <<"dR too small: m="<<m<<" "<<sM<<endl; 
		exit( 0 ); 
	}
	return true;
}	// end checkdR
