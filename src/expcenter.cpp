#include <fstream>
#include <string>
#include "expcenter.h"
#include "lautil.h"
#include "readutil.h"

// create static variables
double CExpCenter::m_sdiel;
double CExpCenter::m_kappa;

double CSolExpCenter::CONST2[N_POLES];
double CSolExpCenter::CONST3[N_POLES];
int CSolExpCenter::id[N_POLES];
int CSolExpCenter::idk[_PSQH_];

#ifdef __ACCURATE__
	const double SOLVE_TOL = ( 1e-8 );
	const int SOLVE_MAX_CT = 100;
#else 
	const double SOLVE_TOL = ( 1e-6 );
	const int SOLVE_MAX_CT = 20;
#endif

/*#########################################################*/
/*#########################################################*/
////////////////////////////////////////////////
// CExpCenter
////////////////////////////////////////////////
/*#########################################################*/
/*#########################################################*/

/******************************************************************/
/******************************************************************/
/**
 * CExpCenter initConstants for initializing constants of a solvent
 exposed expansion.  Also calls an initialization of the RExpan object
 constants
 ******************************************************************/
void 
CExpCenter::initConstants( double sdiel, double kappa )
{
	m_sdiel = sdiel;
	m_kappa = kappa;
	CRExpan::initConstants( kappa );
}	// end initConstants

/*#########################################################*/
/*#########################################################*/
// member functions
/*#########################################################*/
/*#########################################################*/

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

CExpCenter::CExpCenter(  )
		: m_cen( CPnt(0,0,0 )), m_cenRot(CPnt(0,0,0)),  m_rad(0.0), m_idiel(1.0), 
		m_bKappa( false ), m_lscale(1.0),m_mscale(1.0)
{
	m_E = CMulExpan(  );
	m_H = m_E;
}	// end CExpCenter(  )

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

CExpCenter::CExpCenter( int ki, CPnt cen, double rad )
		: m_cen( cen ), m_cenRot(cen), m_rad(rad),m_lscale(rad), m_mscale(1.0/m_rad),
		m_p( N_POLES ), m_id(ki)
{ }

/******************************************************************/
/******************************************************************/
/**
 *  Constructor for ExpCenter class
 ******************************************************************/
CExpCenter::CExpCenter( int ki, CPnt cen, double rad, double idiel, bool bKappa,
		const vector<double> &chgAssigned, const vector<CPnt> &posAssigned, 
		CQuat* orient, CRotCoeff* rot )
		: m_cen( cen ), m_cenRot(cen), m_rad(rad), m_idiel(idiel), m_bKappa(bKappa), 
		m_ch( chgAssigned ), m_pos(posAssigned), m_lscale(rad), 
		m_mscale( 1.0/m_rad ), m_orient(orient), m_rot(rot),
		m_p( N_POLES ), m_id(ki)
{ 
	setFixedCharges(  ); 
} // end CExpCenter(  ..  )

/*#########################################################*/
/*#########################################################*/
// generate m_E, m_H
// m_E, m_rE strictly not neccessary because 
// we depend on m_Efix and m_Lfix of CSolvent centers
// to represent fixed charges that get 
// subsequently incorporated into m_H
// useful only as initial guess for H=E -- to delete this later ???
/*#########################################################*/
/*#########################################################*/

void
CExpCenter::setFixedCharges(  )
{
	// use scaled charges ( values in vacuum )
	vector<REAL> scaledchg;
	double sumcharge = 0.0; // debug
	for( int i=0; i<m_ch.size( ); i++)
	{
		scaledchg.push_back(  m_ch[i]/m_idiel  );
		sumcharge += scaledchg.back(  );
	}
	m_E  = CMulExpan( scaledchg, m_pos, N_POLES, m_bKappa, m_lscale );
	m_rE = m_E;

	//set scale for L
	m_LS.setScale(  m_lscale  );
} // end setFixedCharges

/*#########################################################*/
/*#########################################################*/
//Computes interaction energy of center with ext potential
/*#########################################################*/
/*#########################################################*/

const REAL 
CExpCenter::computePot(  ) const
{
	return  inprod( getH( ),  getLS());
}	// end computePot(  )

/*#########################################################*/
/*#########################################################*/
// force in original reference frame 
// to convert to labframerotate by conj( m_orient ) 
// this will compute force on a specific sphere ( S. Liu )
/*#########################################################*/
/*#########################################################*/

const CPnt
CSolExpCenter::computeForceOn_0(  ) const
{
	if(  IsNoInteractionList( )) return CPnt();

	CPnt force = CPnt(  );

	// spheres that had been polarized
	if(  !IsEmptyInterPolList( ) )
	{
		force  -= inprod(  m_H , m_gLHN ); 
		#ifndef __NOGRADPOL__
			force  -= inprod(  m_gH, m_LS );
		#endif
	}

	// spheres that are not polarized, but have interaction partners
	else
		force  -= inprod(  *m_Hself , m_gLHN );

	if( force.norm( ) > 5 || isnan(force.norm()))
	{
		cout <<"m_H "<<m_H<<endl; 
		cout <<"mgLHN "<<m_gLHN<<endl;
	}

	return force;
}	// end computeForceOn_0

/*#########################################################*/
/*#########################################################*/
// assignment operator overloading ( S. Liu )
/*#########################################################*/
/*#########################################################*/

CExpCenter & 
CExpCenter::operator=( const CExpCenter &M )
{
	m_cen = M.m_cen;
	m_rad = M.m_rad;
	m_idiel = M.m_idiel;
	m_ch = M.m_ch;
	m_pos = M.m_pos;
	m_bKappa = M.m_bKappa;
	m_lscale = M.m_lscale;
	m_mscale = M.m_mscale;

	//m_L = M.m_L;
	m_E = M.m_E;
	m_H = M.m_H;

	return *this;
}	// end CExpCenter::operator=


/*#########################################################*/
////////////////////////////////////////////////
// CSolExpCenter
////////////////////////////////////////////////
/*#########################################################*/

/******************************************************************/
/******************************************************************/
/**
 * CSolExpCenter initConstants for initializing constants of a solvent
 exposed expansion
 ******************************************************************/
void 
CSolExpCenter::initConstants( double sdiel, double kappa )
{
	CExpCenter::initConstants( sdiel, kappa );

	// compute ID, constants
	for( int l=0; l<N_POLES; l++ )
	{
		CONST2[l] = ( 2*l+1 ) / (4.0*M_PI); //The constant (2l+1)/(4*pi) used in EQs 15-17 in Yap 2010
		CONST3[l] = 1 / ( double )(2*l+1);  // 1/(2l+1) 
		id[l] =  l*( l+1 )/ 2;				  // l*(l+1)/2
	}

	// compute IDK: stores index of the first row of every ( l,s ) column
	idk[0] = 0;
	for( int l=1; l<N_POLES; l++ )
	{
		for( int s=0; s <=l; s++ )
		{
			int k = id[l] + s; 
			idk[k] = idk[k-1] + k;
		}
	}
return;
}	// end initConstants

/*#########################################################*/
/*#########################################################*/
// prepare response matrix IMAT
// making use of symmetry in Ynm*Yls
// taken into account conjugation of Ynm with sNeg mNeg
/*#########################################################*/
/*#########################################################*/

void 
CSolExpCenter::build_IMat( const vector<CPnt> &SPE, 
				const vector<CPnt> &SPB, REAL* IMat )
{
	cout <<"static building Imat ... "<<endl;

	int npe = SPE.size(  );
	int npb = SPB.size(  );

	// step 1: compute integrals 
	double *Ylsr_Ynmr = ( double* ) malloc(_PQUADH_*sizeof(double));
	double *Ylsr_Ynmi = ( double* ) malloc(_PQUADH_*sizeof(double));
	double *Ylsi_Ynmr = ( double* ) malloc(_PQUADH_*sizeof(double));
	double *Ylsi_Ynmi = ( double* ) malloc(_PQUADH_*sizeof(double));

	for( int k=0; k < _PQUADH_; k++ )
	{
		Ylsr_Ynmr[k] = 0.0;
		Ylsr_Ynmi[k] = 0.0;
		Ylsi_Ynmr[k] = 0.0;
		Ylsi_Ynmi[k] = 0.0;
	} 

	double start = read_timer(  );

	if( npe <= npb ) 
		computeIntegralE( Ylsr_Ynmr, Ylsr_Ynmi, Ylsi_Ynmr, Ylsi_Ynmi, SPE, SPB );
	else  
		computeIntegralB( Ylsr_Ynmr, Ylsr_Ynmi, Ylsi_Ynmr, Ylsi_Ynmi, SPE, SPB );

	double time1 = read_timer(  );
	//  cout << "Time taken to add integrals "<<time1 - start <<endl;

	// step 2 : populate Imat
	const int len = _PSQ_;
	int k = 0;  // Imat counter ( column major ) 

	// to handle conjugation
	int sNeg = -1;
	int mNeg = -1;

	for( int l=0; l<N_POLES; l++ )
	{
		double ytemp;

		//--------------------------------
		// s=0
		for( int n=0; n<N_POLES; n++ )
		{
			double fact_n = CONST2[n];
			for( int m=0; m<=n; m++ )
			{ 
				bool bUpper = ( n<l || (n==l && m<=0 ) );

				ytemp = (  bUpper? getYY(Ylsr_Ynmr, n,m,l,0 ) : getYY(Ylsr_Ynmr, l,0,n,m)   );
				IMat[k++] = fact_n * ytemp;  // RR-part

				if( m!=0 )
				{
					ytemp = (  bUpper? getYY(Ylsr_Ynmi, n,m,l,0 ) 
									: getYY( Ylsi_Ynmr, l,0,n,m )   );
					IMat[k++] = mNeg * -1.0 * fact_n * ytemp;  //IR-part
				}
			} //m
		}//n

		//--------------------------------
		// s>0
		for( int s=1; s<=l; s++ )
		{
			// ( l,s )real columns	
			for( int n=0; n<N_POLES; n++ )
			{
				double fact_n = CONST2[n];
				for( int m=0; m<=n; m++ )
				{	
					bool bUpper = ( n<l || (n==l && m<=s ) );

					ytemp = (  bUpper? getYY(Ylsr_Ynmr, n,m,l,s ) : 
								getYY( Ylsr_Ynmr, l,s,n,m )   );	  
					IMat[k++] =  2.0 * fact_n * ytemp; //   RR-part	

					if( m!=0 ) 
					{
						ytemp = (  bUpper? getYY(Ylsr_Ynmi, n,m,l,s ) : 
									getYY( Ylsi_Ynmr, l,s,n,m )   );	  
						IMat[k++] = mNeg * -2.0 * fact_n * ytemp; //   IR-part		    
					}
				} //m
			}//n

			// ( l,s )imag columns	
			for( int n=0; n<N_POLES; n++ )
			{
				double fact_n = CONST2[n];
				for( int m=0; m<=n; m++ )
				{
					bool bUpper = ( n<l || (n==l && m<=s ) );

					ytemp = (  bUpper? getYY(Ylsi_Ynmr, n,m,l,s ) : 
									getYY( Ylsr_Ynmi, l,s,n,m )   );	  
					IMat[k++]   = sNeg * -2.0 * fact_n* ytemp; //   RI-part

					if( m!=0 )
					{
						ytemp = (  bUpper? getYY(Ylsi_Ynmi, n,m,l,s ) : 
									getYY( Ylsi_Ynmi, l,s,n,m )   );	  
						IMat[k++] = 2.0 * fact_n * ytemp; //   II-part
					}
				} //m
			}//n
		} //s
	}//l

	double time3 = read_timer(  );cout << "Time taken to prepare Imat "
						<<time3 - start <<endl;

	double vm, rss;
	process_mem_usage( vm, rss );  cout << "After integral: VM: " << vm 
						<< "; RSS: " << rss << endl;

	free( Ylsr_Ynmr );  
	free( Ylsr_Ynmi ); 
	free( Ylsi_Ynmr ); 
	free( Ylsi_Ynmi ); 

	//printMat( IMat );
	//cout <<"ImatE completed ... "<<endl;

	return;
} // end build_IMat

/******************************************************************/
/******************************************************************/
/**
*  for building Imat ( S. Liu )
******************************************************************/
CSolExpCenter::CSolExpCenter( int ki, CPnt cen, double rad,
							 const vector<CPnt> &SPE, const vector<CPnt> &SPB,
							 const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
							 const vector<int> &intraPolList_near,
							 double idiel, bool bKappa, const vector<double> &allchg,
							 const vector<CPnt> &allpos,
							 const vector<double> &chgAssigned, const vector<CPnt> &posAssigned )
: CExpCenter( ki, cen, rad, idiel, bKappa, chgAssigned, posAssigned ),
m_SPE( SPE ),m_SPB(SPB), m_SPx(SPx), m_SPExSize(SPx.size()),
m_nSPx( nSPx ), m_neigh(neigh),
m_intraPolList_near( intraPolList_near ),
m_bRead( false ), m_bOwnMat(true)
{
	printf( "Sphere %d CSolExpCenter by building own IMat\n", m_id );
	initMyConstants(  );
	setFixedCharges( allchg, allpos );
	initSurfaceCharges(  );
	
	printf( "Building Mat for sphere %d\n ",m_id );
	double start = read_timer(  );
	char fname[80];   sprintf( fname, "imat.sp%d.out.bin", m_id );
	m_IMat = new REAL[_PQUAD_];
	build_IMatE(  ); writeMat_Unformatted(m_IMat, N_POLES, fname);
	double end = read_timer(  );
	printf( "Time taken to prepare imat [s] : %f \n", end-start );
	
//	cal_DVFixed(  );
}	// end CSolExpCenter

/*#########################################################*/
/*#########################################################*/
// read in pointer to imat
/*#########################################################*/
/*#########################################################*/
CSolExpCenter::CSolExpCenter( int ki, CPnt cen, double rad,
							 const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
							 const vector<int> &intraPolList_near,
							 double idiel, bool bKappa, const vector<double> &allchg, const vector<CPnt> &allpos,
							 const vector<double> &chgAssigned, const vector<CPnt> &posAssigned, REAL* iMat )
: CExpCenter( ki, cen, rad, idiel, bKappa, chgAssigned, posAssigned ),
m_SPx( SPx ), m_SPExSize(SPx.size()), m_nSPx(nSPx) , m_neigh(neigh),
m_intraPolList_near( intraPolList_near ), m_IMat(iMat), m_bRead(false), m_bOwnMat(false)
{
	initMyConstants(  );
	setFixedCharges( allchg, allpos );
	initSurfaceCharges(  );
	double start = read_timer(  );
	cal_DVFixed(  );
} // end CSolExpCenter

/*#########################################################*/
/*#########################################################*/
// read in pointer to imat, and copy expansions
/*#########################################################*/
/*#########################################################*/
CSolExpCenter::CSolExpCenter( int ki, CPnt cen, double rad,
							 const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
							 const vector<int> &intraPolList_near,
							 double idiel, bool bKappa, const vector<double> &allchg, const vector<CPnt> &allpos,
							 const vector<double> &chgAssigned, const vector<CPnt> &posAssigned, REAL* iMat,
							 const CMulExpan *Fself, const CMulExpan *Hself,
							 const CLocalExpan* LF_intraSelf, const CLocalExpan* LH_intraSelf,
							 const CLocalExpan* LF_intraSelf_far, const CLocalExpan* LH_intraSelf_far,
							 const double *qSolvedFself, const double  *qSolvedHself,
							 const double *totalFself, const double  *totalHself,
							 CQuat* orient, CRotCoeff* rot, bool bGrad )
:
CExpCenter( ki, cen, rad, idiel, bKappa, chgAssigned, posAssigned, orient, rot ),
m_SPx( SPx ), m_SPExSize(SPx.size()), m_nSPx(nSPx), m_neigh(neigh),
m_intraPolList_near( intraPolList_near ),
m_IMat( iMat ),
m_Fself( Fself ), m_LF_intraSelf(LF_intraSelf), m_LH_intraSelf(LH_intraSelf),
m_LF_intraSelf_far( LF_intraSelf_far ), m_LH_intraSelf_far(LH_intraSelf_far),
m_qSolvedFself( qSolvedFself ), m_qSolvedHself(qSolvedHself),
m_totalFself( totalFself ), m_totalHself(totalHself),
m_bRead( true ), m_bOwnMat(false), m_bGrad(bGrad)
{
	assert( getSPx( ).size() == getSPExSize());
	initMyConstants(  );
	setFixedCharges( allchg, allpos );
	m_Hself = Hself;
	m_H = *m_Hself;   m_F = *m_Fself;
	initSurfaceCharges(  );
	cal_DVFixed(  );
}	// end CSolExpCenter

/*#########################################################*/
/*#########################################################*/
// for preparing various moltype quantities
/*#########################################################*/
/*#########################################################*/

CSolExpCenter::CSolExpCenter( int ki, CPnt cen, double rad,
							 const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
							 const vector<int> &intraPolList_near,
							 const CMulExpan *Fself, const CMulExpan *Hself,
							 const double *qSolvedFself, const double  *qSolvedHself,
							 const double *totalFself, const double  *totalHself )
:
CExpCenter( ki, cen, rad ),
m_SPx( SPx ), m_SPExSize(SPx.size()), m_nSPx(nSPx), m_neigh(neigh),
m_intraPolList_near( intraPolList_near ),
m_Fself( Fself ),
m_qSolvedFself( qSolvedFself ), m_qSolvedHself(qSolvedHself), m_totalFself(totalFself),
m_totalHself( totalHself )
{
	m_Hself = Hself;
	m_H = *m_Hself;
	m_F = *m_Fself;
} // end CSolExpCenter

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

CSolExpCenter::CSolExpCenter( int ki, CPnt cen, double rad, const CMulExpan *Hself,
							 const vector<CPnt> &SPx, const int &nSPx, const vector<int> &intraPolList_near )
: CExpCenter( ki, cen, rad ),
m_SPx( SPx ), m_SPExSize(SPx.size()), m_nSPx(nSPx),
m_intraPolList_near( intraPolList_near ),
m_bRead( true ), m_bOwnMat(false)
{
	m_Hself = Hself;
	m_H = *m_Hself;
} // end CSolExpCenter

/******************************************************************/
/******************************************************************/
/**
 *  Initialize constants for building Imat
 ******************************************************************/
void
CSolExpCenter::initMyConstants(  )
{
	vector<double> bessel_K( N_POLES+1, 1.0 ), bessel_I(N_POLES+1, 1.0);
	vector<double> Dcoeff_F( N_POLES ), Vcoeff_F(N_POLES),
	Dcoeff_H( N_POLES ), Vcoeff_H(N_POLES);;
	
	double val = m_kappa * m_rad;					// k*r
	double valsqr = val*val;						// k^2*r^2
	double exp_val = exp( -1.0*val );				// exp(-kr)
	double eps = m_idiel/ m_sdiel;					// e_p / e_s
	CLExpan::BESSEL( bessel_K,val  );				//	khat(k*r)
	CMExpan::BESSEL( bessel_I,val  );				// ihat(k*r)
	
	for( int n=0; n<N_POLES; n++ )
	{
		Dcoeff_H[n] = - exp_val * bessel_K[n];	// -exp(-kr)*khat_n(kr)
		Dcoeff_F[n] = 1.0;
		
		Vcoeff_H[n] =  exp_val * (   n * bessel_K[n]	 // From EQ 22b Yap 2010
								  -(2*n+1 )*bessel_K[n+1] );// -exp(-kr)*(khat_n(kr)*n - (2n+1)*khat_n+1(kr))
		Vcoeff_F[n] =  -eps * n;							 // - e_p/e_s * n
		
		// for F-eqn ( based on Von Neumann )
		m_constH1[n]  = Vcoeff_H[n];                       // H-coeff for F-eqn in 22b Yap 2010
		m_constF1[n]  = Vcoeff_F[n] + ( 2*n+1 );           // F-coeff for F-eqn in 22b Yap 2010
		m_constLF1[n] = -m_rad * eps * n;						// for F-eqn, LF of XF component in 14b Yap 2010
		m_constLH1[n] = m_rad * (  n*bessel_I[n] +			// for F-eqn, LH of XF component in 14b Yap 2010
								 valsqr * bessel_I[n+1] / double(2*n+3 ) );
		
		// for H-eqn ( based on Dirchlet )
		m_constH2[n]  = Dcoeff_H[n] + ( 2*n+1 ) / bessel_I[n]; // H-coeff for H-eqn in 22a Yap 2010
		m_constF2[n]  = Dcoeff_F[n];								 // F-coeff for H-eqn in 22a Yap 2010
		m_constLF2[n] = m_rad;										 // for H-eqn, LF of XH component in 14a Yap 2010
		m_constLH2[n] = - m_rad * bessel_I[n];					 // for H-eqn, LH of XH component in 14a Yap 2010
		
		m_constInvH[n] = bessel_I[n] * ( CONST3[n] );// Inverse of common constant: (ihat_n(kr)/(2n+1))
		m_constChgH[n] = CONST2[n] / bessel_I[n];	 // Common constant in many H equs: (2n+1)/(4pi*ihat_n(kr))
	}
}	/// end initMyConstants


/******************************************************************/
/******************************************************************/
/**
 *  setFixedCharges function. For a given CG sphere, 1. generate m_E, m_H
 2. set charges for current center and generate m_Efix, m_Lfix
 ******************************************************************/
void
CSolExpCenter::setFixedCharges( const vector<double> &allchg,
							   const vector<CPnt> &allpos )
{
	// assigned all fixed charges to inside/outside sphere
	// use scaled charges
	vector<double> chgIn, chgOut;
	vector<CPnt>   posIn, posOut;
	double rsq = m_rad*m_rad;
	
	for( int i=0; i<allpos.size( ); i++)
	{
		if(  allpos[i].normsq( ) < rsq)
		{
			chgIn.push_back(  allchg[i]/m_idiel  ); // q/e_p
			posIn.push_back(  allpos[i]  );
		}
		else
		{
			chgOut.push_back(  allchg[i]/m_idiel  ); // q/e_p
			posOut.push_back(  allpos[i]  );
		}
	}
	
	// generate m_Efix, m_Lfix
	m_Efix = CMulExpan( chgIn, posIn, N_POLES, m_bKappa, m_lscale );
	m_Lfix = CLocalExpan( chgOut, posOut, N_POLES, m_bKappa, m_lscale );
	
	// set scale and range of m_S, LH, LF
	if( !m_bRead )
	{
		m_H = m_rE;
		m_F = CMulExpan( CRange(0, m_p ), m_lscale);
	}
	
	m_Hrot = m_H;
	
	m_LFS_Far = CLocalExpan( CRange(0, m_p ), m_lscale);
	m_LHS_Far = CLocalExpan( CRange(0, m_p ), m_lscale);
	m_LFS = CLocalExpan( CRange(0, m_p ), m_lscale);
	m_LHS = CLocalExpan( CRange(0, m_p ), m_lscale);
}

/*#########################################################*/
/*#########################################################*/
// init after system size ( nmol ) is known
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::initGradient( int nmol )
{
	m_maxgF = 0;
	m_maxgH = 0;
	
	m_gSolvedF.resize( m_SPExSize );
	m_gSolvedH.resize( m_SPExSize );
	
	m_gLHN.reset( m_p );  m_gLHN.setScale(m_lscale);
}	// end initGradient

/******************************************************************/
/******************************************************************/
/**
 *  Function to initialize surface charges for self-polarization
 ******************************************************************/
void
CSolExpCenter::initSurfaceCharges(  )
{
	m_qSolvedF.resize( m_SPExSize );
	m_qSolvedH.resize( m_SPExSize );
	
	m_dA = 4*M_PI/double(  m_nSPx  );
	
	computeExposedSurfaceChargesH(  ); // either read-in or initialized to rE
	if( m_bRead )
		computeExposedSurfaceChargesF(  );
	else
		m_totalF = 0;
}	// end initSurfaceCharges

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

double
CSolExpCenter::solveSurfaceCharges(  )
{
	solveSurfaceCharges( m_F, m_H, false );
}	// end solveSurfaceCharges

/*#########################################################*/
/*#########################################################*/
// calculates the surface charges iteratively using y=Ax
// based on molecule in initial orientation
//( i.e. LH must be rotated by to
//initial position by -rot, and keep
// a copy of H ( as Hrot ) to polarize other spheres
/*#########################################################*/
/*#########################################################*/

double
CSolExpCenter::solveSurfaceCharges( const CMulExpan &Finit,
								   const CMulExpan &Hinit, bool bRotateH )
{
	double devH;
	const int len = m_p * m_p;
	
	vector<double> F = Finit.getVector( m_p );
	vector<double> H = Hinit.getVector( m_p );
	
	devH = solveFH( F, H, m_LFS.getVector(m_p ), m_LHS.getVector(m_p), m_p, m_p, 1);
	
	m_F = CMulExpan( F, CRange(0, m_p ), m_mscale);
	m_H = CMulExpan( H, CRange(0, m_p ), m_mscale);
	if( bRotateH )
		m_rot->rotateWithXi( m_H,m_Hrot , 1, m_p, true ); // rotate for mutual pol; not self pol
	
	//   cout <<"F"<<m_F<<endl;
	// cout <<"H"<<m_H<<endl;
	m_dev = devH;
	
	// also recalculate qsolved ( exposed surface charges )
	computeExposedSurfaceChargesFH(  );
	
	return devH;
}	// end solveSurfaceCharges

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

double
CSolExpCenter::solveSurfaceGradient( int j, CGradExpan &gF, CGradExpan &gH, CGradExpan &gHrot,
									const CGradExpan & LGF_j, const CGradExpan & LGH_j )
{
	const int len = m_p * m_p;
	
	vector<double> gF_j = gF.getOneVector(  );
	vector<double> gH_j = gH.getOneVector(  );
	
	double devH = solveFH( gF_j, gH_j, LGF_j.getOneVector( ), LGH_j.getOneVector(), m_p, m_p, 3);
	
	gF.setFromOneVector( gF_j );
	gH.setFromOneVector( gH_j );
	m_rot->rotateWithXi( gH, gHrot , 1, m_p, true );
	m_dev = devH;
	
	computeExposedSurfaceGradientFH( gF, gH );
	
	return devH;
}	// end solveSurfaceGradient

/******************************************************************/
/******************************************************************/
/**
 *  compute exposed surface integral directly using EXPOSED surface points
 ******************************************************************/
void
CSolExpCenter::computeIntegralE( double * Yrr,double * Yri,double * Yir,double * Yii ,
								const vector<CPnt> & SPE, const vector<CPnt> & SPB )
{
	cout <<"compute integral using exposed points"<<endl;
	
	int npe = SPE.size(  );
	int npb = SPB.size(  );
	int np = npb + npe;
	
	for( int h=0; h<npe; h++ )
	{
		double w=0;
		// convert the position relative to the center to spherical coordinates.
		CSpPnt q = CartToSph( SPE[h] );
		CSHExpan Y( q.theta( ), q.phi(), N_POLES);
		
		// collect sums for ( n,m ) rows x (l,s) column
		int kY=0;
		for( int l=0; l<N_POLES; l++ )
		{
			for( int s=0; s<=l; s++ )
			{
				double Ylsr, Ylsi;
				if( s==0 )
				{
					Ylsr = 	Y( l,0 );
					Ylsi =  0.0;
				}
				else
				{
					Ylsr = Y( l, 2*s-1 );
					Ylsi = Y( l, 2*s ) ;
				}				
				for( int n=0; n<=l; n++ )
				{
					for( int m=0; m<=n; m++ )
					{
						if( n==l  && m > s ) break;
						
						double Ynmr, Ynmi;
						if( m==0 )
						{
							Ynmr = Y( n,0 );
							Ynmi = 0.0;
						}
						else
						{
							Ynmr = Y( n, 2*m-1 ) ;
							Ynmi = Y( n, 2*m ) ;
						}
						
						// integrate using the appropriate integration rules
						if( m==0 && s==0 && (n+l )%2==0) //simpson's rule
						{
							if( h > 0 && h < np )
							{
								if( h % 2 == 0  )
									w = 2.0/3.0; // even
								else w = 4.0/3.0; //odd
							}
							else w = 1.0/3.0;
						}
						else // rectangular
							w = 1.0;
						
						Yrr[kY] += w * Ylsr*Ynmr;
						Yri[kY] += w * Ylsr*Ynmi;
						Yir[kY] += w * Ylsi*Ynmr;
						Yii[kY] += w * Ylsi*Ynmi;
						
						kY++;
					}// m
				}// n
			}// s
		}// l
	}//h
	
	// clean up integral, using the fact that when m=s the integral is real
	// divide the sums by to get surface integral
	double dA = 4*M_PI / ( double )(np-1);
	
	for( int l=0, kY=0; l<N_POLES; l++ )
	{
		for( int s=0; s<=l; s++ )
		{
			for( int n=0; n<=l; n++ )
			{
				for( int m=0; m<=n; m++ )
				{
					if(  n==l && m > s  ) break;
					
					if(  n==l && m==s  )
					{
						if( m==0 )
						{
							Yrr[kY] *= dA;
							Yii[kY] = 0.0;
						}
						else
						{
							Yrr[kY] *= dA;
							Yii[kY] *= dA;
						}
						Yri[kY]  *= dA;
						Yir[kY]  *= dA;
					}
					else
					{
						Yrr[kY] *= dA;
						Yri[kY] *= dA;
						Yir[kY] *= dA;
						Yii[kY] *= dA;
					}
					kY++;
				}
			}
		}
	}
}	// end computeIntegralE


/******************************************************************/
/******************************************************************/
/**
 * compute exposed surface integral indirectly using BURIED surface points
 ******************************************************************/
void
CSolExpCenter::computeIntegralB( double * Yrr,double * Yri,double * Yir,double * Yii,
								const vector<CPnt> & SPE, const vector<CPnt> & SPB  )
{
	cout <<"compute integral using buried points"<<endl;
	int npe = SPE.size(  );
	int npb = SPB.size(  );
	int np = npb + npe;
	
	for( int h=0; h<npb; h++ )
	{
		double w=0;
		// convert the position relative to the center to spherical coordinates.
		CSpPnt q = CartToSph( SPB[h] );
		CSHExpan Y( q.theta( ), q.phi(), N_POLES);
		
		// collect sums for ( n,m ) rows x (l,s) column
		for( int l=0, kY=0; l<N_POLES; l++ )
		{
			for( int s=0; s<=l; s++ )
			{
				double Ylsr, Ylsi;
				if( s==0 )
				{
					Ylsr = 	Y( l,0 );
					Ylsi =  0.0;
				}
				else
				{
					Ylsr = Y( l, 2*s-1 );
					Ylsi = Y( l, 2*s ) ;
				}
				
				for( int n=0; n<=l; n++ )
				{
					for( int m=0; m<=n; m++ )
					{
						if( n==l  && m > s ) break;
						
						double Ynmr, Ynmi;
						if( m==0 )
						{
							Ynmr = Y( n,0 );
							Ynmi = 0.0;
						}
						else
						{
							Ynmr = Y( n, 2*m-1 ) ;
							Ynmi = Y( n, 2*m ) ;
						}
						
						// integrate using the appropriate integration rules
						if( m==0 && s==0 && (n+l )%2==0) //simpson's rule
						{
							if( h > 0 && h < np )
							{
								if( h % 2 == 0  ) w = 2.0/3.0; // even
								else w = 4.0/3.0; //odd
							}
							else w = 1.0/3.0;
						}
						else // rectangular
							w = 1.0;
						
						Yrr[kY] += w * Ylsr*Ynmr;
						Yri[kY] += w * Ylsr*Ynmi;
						Yir[kY] += w * Ylsi*Ynmr;
						Yii[kY] += w * Ylsi*Ynmi;
						
						kY++;
					}//nm
				}
			}// s
		}// l
		
		if(  h % 10000 == 0  ) cout <<"completed "<<h<<" points "<<endl;
	}//h
	
	// modify to get exposed integral
	// divide the sums by to get surface integral
	double dA = 4*M_PI / ( double )(np-1);
	
	for( int l=0, kY=0; l<N_POLES; l++ )
	{
		double delta = 1.0/CONST2[l];
		for( int s=0; s<=l; s++ )
		{
			for( int n=0; n<=l; n++ )
			{
				for( int m=0; m<=n; m++ )
				{
					if(  n==l && m > s  ) break;
					
					if(  n==l && m==s  )
					{
						if( m==0 )
						{
							Yrr[kY] = delta - dA*Yrr[kY];
							Yii[kY] = 0.0;
						}
						else
						{
							Yrr[kY] = 0.5*delta - dA*Yrr[kY];
							Yii[kY] = 0.5*delta - dA*Yii[kY];
						}
						Yri[kY] *= -dA;
						Yir[kY] *= -dA;
					}
					else
					{
						Yrr[kY] *= -dA;
						Yri[kY] *= -dA;
						Yir[kY] *= -dA;
						Yii[kY] *= -dA;
					}
					kY++;
				}
			}
		}
	}
} // end computeIntegralB

/******************************************************************/
/******************************************************************/
/**
 *  Function to prepare response matrix IMatE making use of
 symmetry in Ynm*Yls taken into account conjugation of Ynm with sNeg mNeg
 ******************************************************************/
void
CSolExpCenter::build_IMatE(  )
{
	cout <<"building ImatE  ... "<<endl;
	
	const int npe = m_SPE.size(  );
	const int npb = m_SPB.size(  );
	
	// step 1: compute integrals
	double *Ylsr_Ynmr = ( double* ) malloc(_PQUADH_*sizeof(double));
	double *Ylsr_Ynmi = ( double* ) malloc(_PQUADH_*sizeof(double));
	double *Ylsi_Ynmr = ( double* ) malloc(_PQUADH_*sizeof(double));
	double *Ylsi_Ynmi = ( double* ) malloc(_PQUADH_*sizeof(double));
	
	for( int k=0; k < _PQUADH_; k++ )
	{
		Ylsr_Ynmr[k] = 0.0;
		Ylsr_Ynmi[k] = 0.0;
		Ylsi_Ynmr[k] = 0.0;
		Ylsi_Ynmi[k] = 0.0;
	}
	
	double start = read_timer(  );
	
	if( npe <= npb ) // If there are less exposed surface points than buried
		computeIntegralE( Ylsr_Ynmr, Ylsr_Ynmi, Ylsi_Ynmr, // directly compute their integrals
						 Ylsi_Ynmi, m_SPE, m_SPB );
	else
		computeIntegralB( Ylsr_Ynmr, Ylsr_Ynmi, Ylsi_Ynmr,
						 Ylsi_Ynmi, m_SPE, m_SPB );
	
	double time1 = read_timer(  );
	cout << "Time taken to add integrals "<<time1 - start <<endl;
	
	// step 2 : populate Imat
	const int len = _PSQ_;
	vector<REAL> fact( N_POLES );
	
	int k = 0;  // Imat counter ( column major )
	
	// to handle conjugation
	int sNeg = -1;
	int mNeg = -1;
	
	for( int l=0; l<N_POLES; l++ )
	{
		double ytemp;		
		//--------------------------------
		// s=0
		for( int n=0; n<N_POLES; n++ )
		{
			double fact_n = CONST2[n];
			for( int m=0; m<=n; m++ )
			{
				bool bUpper = ( n<l || (n==l && m<=0 ) );
				ytemp = (  bUpper? getYY(Ylsr_Ynmr, n,m,l,0 ) : getYY(Ylsr_Ynmr, l,0,n,m)   );
				m_IMat[k++] = fact_n * ytemp;  // RR-part
				
				if( m!=0 )
				{
					ytemp = (  bUpper? getYY(Ylsr_Ynmi, n,m,l,0 ) : getYY(Ylsi_Ynmr, l,0,n,m)   );
					m_IMat[k++] = mNeg * -1.0 * fact_n * ytemp;  //IR-part
				}
			} //m
		}//n
		
		//--------------------------------
		// s>0
		for( int s=1; s<=l; s++ )
		{
			// ( l,s )real columns
			for( int n=0; n<N_POLES; n++ )
			{
				double fact_n = CONST2[n];
				for( int m=0; m<=n; m++ )
				{
					bool bUpper = ( n<l || (n==l && m<=s ) );
					ytemp = (  bUpper? getYY(Ylsr_Ynmr, n,m,l,s ) : getYY(Ylsr_Ynmr, l,s,n,m)   );
					m_IMat[k++] =  2.0 * fact_n * ytemp; //   RR-part
					
					if( m!=0 )
					{
						ytemp = (  bUpper? getYY(Ylsr_Ynmi, n,m,l,s ) : getYY(Ylsi_Ynmr, l,s,n,m)   );
						m_IMat[k++] = mNeg * -2.0 * fact_n * ytemp; //   IR-part
					}
				} //m
			}//n
			
			// ( l,s )imag columns
			for( int n=0; n<N_POLES; n++ )
			{
				double fact_n = CONST2[n];
				for( int m=0; m<=n; m++ )
				{
					bool bUpper = ( n<l || (n==l && m<=s ) );
					ytemp = (  bUpper? getYY(Ylsi_Ynmr, n,m,l,s ) : getYY(Ylsr_Ynmi, l,s,n,m)   );
					m_IMat[k++]   = sNeg * -2.0 * fact_n* ytemp; //   RI-part
					
					if( m!=0 )
					{
						ytemp = (  bUpper? getYY(Ylsi_Ynmi, n,m,l,s ) : getYY(Ylsi_Ynmi, l,s,n,m)   );
						m_IMat[k++] = 2.0 * fact_n * ytemp; //   II-part
					}
				} //m
			}//n
			
		} //s
	}//l
	
	double time3 = read_timer(  );
	cout << "Time taken to prepare Imat "<<time3 - start <<endl;
	
	double vm, rss;
	process_mem_usage( vm, rss );
	cout << "After integral: VM: " << vm << "; RSS: " << rss << endl;
	
	free( Ylsr_Ynmr );
	free( Ylsr_Ynmi );
	free( Ylsi_Ynmr );
	free( Ylsi_Ynmi );
	return;
} // end build_IMatE


/******************************************************************/
/******************************************************************/
/**
 *  Function to caclulate the fixed coefficients for dirichelet and
 van neumann solutions. Used in EQ 14 of Yap 2010
 ******************************************************************/
void
CSolExpCenter::cal_DVFixed(  )
{
	const double eps = m_idiel/m_sdiel;
	for( int n=0, k=0; n<N_POLES; n++ )
	{
		const double Dcoeff_E =  1.0;					// 1
		const double Dcoeff_L =  m_rad;				// rho
		const double Vcoeff_E =  eps * ( n+1 );		//	e_p/e_s * (n+1)
		const double Vcoeff_L = -eps * n * m_rad;	// -e_p/e_s * n * rho
		
		for( int mm=0; mm<2*n+1; mm++, k++ )
		{
			m_Dfix[k] = Dcoeff_E * m_Efix[k] + Dcoeff_L * m_Lfix[k]; // Enm + rho*LE
			m_Vfix[k] = Vcoeff_E * m_Efix[k] + Vcoeff_L * m_Lfix[k]; // e_p/e_s*(n+1)*Enm - e_p/e_s*n*rho*LE
		}
	}
	return;
}	// end cal_DVFixed

/*#########################################################*/
/*#########################################################*/
// iteratively solves for F, H using D, V separately
// conjugation taken into account in IMAT
// I think this is the core function for doing mutual polarization
// In general, what function 'solveFH' did
// is to calculate F and H
// iteratively with given input of F, H, LF and LH.
// It has called function applyMM, which can
// be found in file 'lautil.cpp',
// to do matrix vector products.
// So the job 'solveFH' has done is to
// solve linear equations in 4a and 4b
// in 2013 JCTC paper of THG and Enghui Yap.
//
// Note that 'solveFH' is not only used
// for solving multiple expansion of
// surface charge density like F and H,
// but also the gradients of multiple
// Add a comment to this line
// expansions, such as equation 11a and 11b in 2013
// JCTC paper of THG and Enghui Yap. ( S. Liu )
/*#########################################################*/
/*#########################################################*/

double
CSolExpCenter::solveFH( vector<double> &F, vector<double>  &H,
					   const vector<double> &LF, const vector<double>  &LH,
					   int pm, int pn, int D ) const
{
	assert( D==1 || D==3 );
	const int p2 = pm*pm;
	const int mlen = D*pm*pm; // nm
	const int nlen = D*pn*pn; // ls
	
	const int expectedLength = p2*D;
	assert(  F.size( ) == expectedLength);
	assert(  H.size( ) == expectedLength);
	
	vector<double> hbase( mlen ), fbase(mlen);
	vector<double> hx( nlen ), fx(nlen);
	vector<double> oldH( mlen ),oldF(mlen), outerOldH(mlen),outerOldF(mlen);
	
	// compute hbase and fbase to incorporate LF, LH
	if( D==1 )
		compute_FHbase( &(LH[0] ), &(LF[0]), &(fbase[0]), &(hbase[0]), pm);
	else
		compute_FHbase3( &(LH[0] ), &(LF[0]), &(fbase[0]), &(hbase[0]), pm);
	
	// initialize old value ( just H since we choose to use it for convergence )
	oldH = H; outerOldH = H;
	
	// start polarizing
	double dev, dev_h, dev_f;
	
	// 1st cycle
	compute_fx( &(fx[0] ), &(fbase[0]), &(F[0]), &(H[0]), pn, D);
	applyMMat( m_IMat, &(fx[0] ), &(F[0]), 1.0, 0.0, p2, D, p2);
	revertF( &(F[0] ), pm, D);
	
	compute_hx( &(hx[0] ), &(hbase[0]), &(F[0]), &(H[0]), pn, D);
	applyMMat( m_IMat, &(hx[0] ), &(H[0]), 1.0, 0.0, p2, D, p2);
	revertH( &(H[0] ), pm, D);
	dev_h = computeDev( &(oldH[0] ), &(H[0]), pm,D); oldH = H;
	
	dev = dev_h;
	
	int ct = 1;
	while( dev > SOLVE_TOL && ct < SOLVE_MAX_CT )
	{
		compute_fx( &(fx[0] ), &(fbase[0]), &(F[0]), &(H[0]), pn, D);
		applyMMat( m_IMat, &(fx[0] ), &(F[0]), 1.0, 0.0, p2, D, p2);
		revertF( &(F[0] ), pm, D);
		
		compute_hx( &(hx[0] ), &(hbase[0]), &(F[0]), &(H[0]), pn, D);
		applyMMat( m_IMat, &(hx[0] ), &(H[0]), 1.0, 0.0, p2, D, p2);
		revertH( &(H[0] ), pm, D);
		dev_h = computeDev( &(oldH[0] ), &(H[0]), pm,D); oldH = H;
		
		dev = dev_h;
		ct++;
	}
	
	// calculate deviation from old m_H, set new m_F, m_H
	double devH = computeDev2(  &(outerOldH[0] ), &(H[0]), pm, D);	
	return devH;
} // end solveFH

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::revertF( double * F, int p, int D ) const
{
	for( int d=0, k=0; d<D; d++ )
	{
		for( int n=0; n<p; n++ )
		{
			const double fact = ( CONST3[n] );
			for( int mm=0; mm<2*n+1; mm++, k++ )
				F[k] *= fact;
		}
	}
}	// end revertF

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::revertH( double * H, int p, int D ) const
{
	for( int d=0, k=0; d<D; d++ )
	{
		for( int n=0; n<p; n++ )
		{
			const double fact = m_constInvH[n];
			for( int mm=0; mm<2*n+1; mm++, k++ )
				H[k] *= fact;
		}
	}
} // end revertH

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::compute_FHbase( const double* LH, const double *LF,
							  double* fbase,double* hbase, int pm ) const
{
	for( int n=0, k=0; n<pm; n++ )
	{
		const double cLH1 = m_constLH1[n];
		const double cLH2 = m_constLH2[n];
		const double cLF1 = m_constLF1[n];
		const double cLF2 = m_constLF2[n];
		
		for( int mm=0; mm<2*n+1; mm++, k++ )
		{
			fbase[k] = m_Vfix[k] + cLH1 * LH[k] + cLF1 * LF[k];
			hbase[k] = m_Dfix[k] + cLH2 * LH[k] + cLF2 * LF[k];
		}
	}
	
	return;
}	// end compute_FHbase

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::compute_fx( double *fx, const double *fbase, const double *F,
						  const double *H, int p, int D ) const
{
	for( int d=0, k=0; d<D; d++ )
	{
		for( int n=0; n<p; n++ )
		{
			const double cH1 = m_constH1[n];
			const double cF1 = m_constF1[n];
			for( int mm=0; mm<2*n+1; mm++, k++ )
				fx[k] = fbase[k] +  cH1 * H[k] + cF1 * F[k];
		}
	}
}	// end compute_fx


/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::compute_hx( double *hx, const double *hbase, const double *F,
						  const double *H, int p, int D ) const
{
	for( int d=0, k=0; d<D; d++ )
	{
		for( int n=0; n<p; n++ )
		{
			const double cH2 = m_constH2[n];
			const double cF2 = m_constF2[n];
			for( int mm=0; mm<2*n+1; mm++, k++ )
				hx[k] = hbase[k] +  cH2 * H[k] + cF2 * F[k];
		}
	}
}	// end compute_hx

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::compute_FHbase3( const double* LH, const double *LF,
							   double* fbase,double* hbase, int pm ) const
{
	for( int d=0, k=0; d<3; d++ )
	{
		for( int n=0; n<pm; n++ )
		{
			const double cLH1 = m_constLH1[n];
			const double cLH2 = m_constLH2[n];
			const double cLF1 = m_constLF1[n];
			const double cLF2 = m_constLF2[n];
			
			for( int mm=0; mm<2*n+1; mm++, k++ )
			{
				fbase[k] = cLH1 * LH[k] + cLF1 * LF[k];
				hbase[k] = cLH2 * LH[k] + cLF2 * LF[k];
			}
		}
	}
	
	return;
}	// end compute_FHbase3


/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

REAL
CSolExpCenter::computeDev( const double *M1, const double *M2, int p, int D )
{
	assert( p>0 );
	REAL sum = 0;
	Complex tM1,tM2;
	
	int k = 0;
	for ( int d = 0; d < D; d++ )
	{
		for ( int n = 0; n < p; n++ )
		{
			for ( int m = 0; m <= n; m++ )
			{
				if( m==0 )
				{
					tM1 = Complex( M1[k],0 );
					tM2 = Complex( M2[k],0 );
					k++;
				}
				else
				{
					tM1 = Complex( M1[k], M1[k+1] );
					tM2 = Complex( M2[k], M2[k+1] );
					k+=2;
				}
				
				REAL s;
				
				if ( fabs(tM1.real( )) < PREC_LIMIT &&
					fabs( tM1.imag( )) < PREC_LIMIT)
				{
					if ( fabs(tM2.real( )) < PREC_LIMIT &&
						fabs( tM2.imag( )) < PREC_LIMIT)
						continue;
					else
						s = 1.0;
				}
				else
				{
					if ( fabs(tM2.real( )) < PREC_LIMIT &&
						fabs( tM2.imag( )) < PREC_LIMIT)
						s = 1.0;
					else
					{
						Complex diff = ( tM1 - tM2 );
						REAL top = diff.real(  )*diff.real() +  diff.imag()*diff.imag();
						REAL bot = tM1.real(  )*tM1.real() +  tM1.imag()*tM1.imag()
						+ tM2.real(  )*tM2.real() +  tM2.imag()*tM2.imag();
						
						if( m!=0 ) s = 2* top / bot; // for m!=0 to account for m<0 conjugates
						else s = top / bot;
						//cout <<"case 4 : s = "<<s<<endl;
					}
				}
				sum += s;
			} // m
		} // n
	} // d
	
	sum /= 4*p*p;
	return sum;
}	// end computeDev

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

REAL
CSolExpCenter::computeDev2( const double *M1, const double *M2, int p, int D )
{
	assert( p>0 );
	
	REAL totdiff=0, totsum = 0;
	Complex tM1,tM2, diff;
	
	for ( int d = 0, k = 0; d < D; d++ )
	{
		for ( int n = 0; n < p; n++ )
		{
			tM1 = Complex( M1[k],0 );
			tM2 = Complex( M2[k],0 );
			k++;
			
			diff = ( tM1 - tM2 );
			totdiff  += diff.real(  )*diff.real() +  diff.imag()*diff.imag();
			totsum   += tM1.real(  )*tM1.real() +  tM1.imag()*tM1.imag()
			+ tM2.real(  )*tM2.real() +  tM2.imag()*tM2.imag();
			
			for ( int m = 1; m <= n; m++ )
			{
				tM1 = Complex( M1[k], M1[k+1] );
				tM2 = Complex( M2[k], M2[k+1] );
				k+=2;
				
				diff = ( tM1 - tM2 );
				totdiff += 2*( diff.real( )*diff.real() +  diff.imag()*diff.imag());
				totsum  += 2*( tM1.real( )*tM1.real() +  tM1.imag()*tM1.imag()
							  + tM2.real(  )*tM2.real() +  tM2.imag()*tM2.imag());
			}// end m
		}// end n
	}// end d
	
	return totdiff / ( 4*totsum ) / double(D);
} // end computeDev2

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::computeExposedSurfaceGradientFH( const CGradExpan &gF,
											   const CGradExpan &gH )
{
	REAL total_gFx = 0.0,total_gFy = 0.0, total_gFz = 0.0;
	REAL total_gHx = 0.0,total_gHy = 0.0, total_gHz = 0.0;
	
	for( int h=0; h<m_SPExSize; h++ )
	{
		CSpPnt q = CartToSph( m_SPx[h] );
		CSHExpan Y( q.theta( ), q.phi(), N_POLES);
		
		CExpan YI;
		
		// F
		YI = Y * CONST2;
		REAL gFx = inprod_unitScale( gF[0], YI ) * m_dA ;
		REAL gFy = inprod_unitScale( gF[1], YI ) * m_dA ;
		REAL gFz = inprod_unitScale( gF[2], YI ) * m_dA ;
		m_gSolvedF[h] = CPnt( gFx, gFy, gFz );
		total_gFx += fabs( gFx ); total_gFy += fabs(gFy); total_gFz += fabs(gFz);
		
		// H
		YI = Y * m_constChgH;
		REAL gHx = inprod_unitScale( gH[0], YI ) * m_dA ;
		REAL gHy = inprod_unitScale( gH[1], YI ) * m_dA ;
		REAL gHz = inprod_unitScale( gH[2], YI ) * m_dA ;
		m_gSolvedH[h] = CPnt( gHx, gHy, gHz );
		total_gHx += fabs( gHx ); total_gHy += fabs(gHy); total_gHz += fabs(gHz);
	}
	
	m_maxgF = ( total_gFx > total_gFy ? total_gFx : total_gFy );
	if( total_gFz > m_maxgF )
		m_maxgF = total_gFz;
	m_maxgH = ( total_gHx > total_gHy ? total_gHx : total_gHy );
	if( total_gHz > m_maxgH )
		m_maxgH = total_gHz;
} // end computeExposedSurfaceGradientFH

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::computeExposedSurfaceChargesFH(  )
{
	m_totalF = 0.0;
	m_totalH = 0.0;
	
	for( int h=0; h<m_SPExSize; h++ )
	{
		CSpPnt q = CartToSph( m_SPx[h] );
		CSHExpan Y( q.theta( ), q.phi(), N_POLES);
		
		CExpan YI;
		
		YI = Y * CONST2;
		m_qSolvedF[h] = inprod_unitScale( m_F, YI ) * m_dA ; // charge = surface charge * area
		m_totalF += fabs( m_qSolvedF[h] );
		
		YI = Y * m_constChgH;
		m_qSolvedH[h] = inprod_unitScale( m_H, YI ) * m_dA ; // charge = surface charge * area
		m_totalH += fabs( m_qSolvedH[h] );
	}
	
	//  assert( fabs(m_totalH ) < 1e6 );
}	// end computeExposedSurfaceChargesFH

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::computeExposedSurfaceChargesFH_( const vector<CPnt> &SPx,
											   const CMulExpan & F, const CMulExpan & H,
											   double *constH,
											   double nSPx,
											   vector<double> &qSolvedF, vector<double> &qSolvedH,
											   double &totalF, double &totalH, int p )
{
	totalF = 0.0;
	totalH = 0.0;
	
	int SPExSize = SPx.size(  );
	qSolvedF.resize( SPExSize );
	qSolvedH.resize( SPExSize );
	
	REAL dA = 4.0*M_PI/double(  nSPx  );
	for( int h=0; h<SPExSize; h++ )
	{
		CSpPnt q = CartToSph( SPx[h] );
		CSHExpan Y( q.theta( ), q.phi(), p);
		
		CExpan YI;
		
		YI = Y * CONST2;
		qSolvedF[h] = inprod_unitScale( F, YI ) * dA ; // charge = surface charge * area
		totalF += fabs( qSolvedF[h] );
		
		YI = Y * constH;
		qSolvedH[h] = inprod_unitScale( H, YI ) * dA ; // charge = surface charge * area
		totalH += fabs( qSolvedH[h] );
	}
	
	//  assert( fabs(m_totalH ) < 1e6 );
}	// end computeExposedSurfaceChargesFH_

/******************************************************************/
/******************************************************************/
/**
 *  Function to initialize H surface charges for self-polarization
 ******************************************************************/
void
CSolExpCenter::computeExposedSurfaceChargesH(  )
{
	m_totalH = 0.0;
	
	for( int h=0; h<m_SPExSize; h++ )
	{
		CSpPnt q = CartToSph( m_SPx[h] );   // Generate spherical coord for surface point
		CSHExpan Y( q.theta( ), q.phi(), N_POLES);  // Generate spherical harmonic expansion Y of EQ 10a Yap 2010
		CExpan YI = Y * m_constChgH;		 // YI = Y*(2n+1)/(4pi*ihat_n(kr))
		
		m_qSolvedH[h] = inprod_unitScale( m_H, YI ) * m_dA ; // charge = surface charge * area
		
		m_totalH += fabs( m_qSolvedH[h] );	// sum of all h = sum(m_H*YI x dA)
	}
}	// end computeExposedSurfaceChargesH

/******************************************************************/
/******************************************************************/
/**
 *  Function to initialize F surface charges for self-polarization
 ******************************************************************/
void
CSolExpCenter::computeExposedSurfaceChargesF(  )
{
	m_totalF = 0.0;
	
	for( int h=0; h<m_SPExSize; h++ )
	{
		CSpPnt q = CartToSph( m_SPx[h] );
		CSHExpan Y( q.theta( ), q.phi(), N_POLES);
		CExpan YI = Y * CONST2;
		
		m_qSolvedF[h] = inprod_unitScale( m_F, YI ) * m_dA ; // charge = surface charge * area
		
		m_totalF += fabs( m_qSolvedF[h] );
	}

}	// end computeExposedSurfaceChargesF

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::getSurfaceChargesH( const vector<CPnt> &SP, vector<double> &qS )
{
	int np = SP.size(  );
	qS.resize( np );
	
	for( int h=0; h<np; h++ )
	{
		CSpPnt q = CartToSph( SP[h] );
		CSHExpan Y( q.theta( ), q.phi(), N_POLES);
		CExpan YI = Y * m_constChgH;
		
		qS[h] = inprod_unitScale( m_H, YI ) * m_dA ; // charge = surface charge * area
	}
}	// getSurfaceChargesH

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::getSurfaceChargesF( const vector<CPnt> &SP, vector<double> &qS )
{
	int np = SP.size(  );
	qS.resize( np );
	
	for( int h=0; h<np; h++ )
	{
		CSpPnt q = CartToSph( SP[h] );
		CSHExpan Y( q.theta( ), q.phi(), N_POLES);
		CExpan YI = Y * CONST2;
		
		qS[h] = inprod_unitScale( m_F, YI ) * m_dA ; // charge = surface charge * area
	}
}	// emd getSurfaceChargesF

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

const CPnt
CSolExpCenter::computeTorqueOn_0( int i ) const
{
	if(  IsEmptyInterPolList( )) return CPnt();
	
	CPnt torque = CPnt(  );
	
	for( int h=0; h<m_SPExSize; h++ )
	{
		CSpPnt q = CartToSph( m_SPx[h] );
		CSHExpan Y( q.theta( ), q.phi(), m_p);
		CExpan YI = Y * m_constChgH;
		
		REAL Yh_LHN = inprod_unitScale( YI, m_LS );
		CPnt Yh_gLHN = inprod_unitScale( YI, m_gLHN );
		
		torque -= cross( m_SPx[h], ( Yh_LHN * m_gSolvedH[h] + m_qSolvedH[h] * Yh_gLHN  )  );
	}
	
	return m_dA * torque;
}	// end computeTorqueOn_0

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

bool
CSolExpCenter::isOnInterPolList( int j, int kj ) const
{
	for( int k=0; k<m_interPolList.size( ); k++)
	{
		if(  m_interPolList[k].mid( ) == j && m_interPolList[k].kid() == kj )
			return true;
	}
	return false;
}

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::rotateCenters(  )
{
	// rotate the position of ki's rotated center wrt rcen
	m_cenRot = *m_orient * m_cen;
	return;
}	// end rotateCenters

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::rotateHself(  )
{
	// for reexpand external
	m_rot->rotateWithXi( *m_Hself, m_HselfRot, 1, m_p, true ); // for recomputeFromSelfVal
	m_Hrot = m_HselfRot; // for recompute, in case Hrot is not ready in parallel runs
	
	return;
}	// end rotateHself

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

#ifdef __NOPOL_UNCHANGED__
void
CSolExpCenter::rotateCurrentH(  )
{
	// for reexpand external
	m_rot->rotateWithXi( m_H, m_Hrot, 1, m_p, true );
	
	return;
}

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

CGradExpan
CSolExpCenter::rotateGH(  ) const
{
	CGradExpan gHrot;
	m_rot->rotateWithXi( m_gH, gHrot , 1, m_p, true );
	return gHrot;
}

#endif

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

CSolExpCenter&
CSolExpCenter::operator=( const CSolExpCenter &M )
{
	( *this ).CExpCenter::operator=(M);
	
	m_SPE = M.m_SPE;
	m_Lfix = M.m_Lfix;
	m_Efix = M.m_Efix;
	//m_LH = M.m_LH;
	//m_LF = M.m_LF;
	m_F = M.m_F;
}	//end CSolExpCenter::operator=

/*#########################################################*/
/*#########################################################*/
//====================
// utility functions
//====================
/*#########################################################*/
/*#########################################################*/

//#ifdef __DEBUG
void
CSolExpCenter::printY( double *Y )
{
	double out;
	for( int n=0; n<N_POLES; n++ )
	{
		for( int m=0; m<=n; m++ )
		{
			for( int l=0; l<N_POLES; l++ )
			{
				for( int s=0; s<=l; s++ )
				{
					if( n > l || (n==l  && m > s )) out = 0.0;
					else
					{
						out = Y[ IDK( l,s ) + id[n] + m];
						if( fabs(out )<1e-5) out = 0.0;
					}
					cout <<out<<" ";
				} //l
			} //s
			cout<<endl;
		}// m
	}// n
}	// end printY

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::printYFull( double *Y )
{
	double out;
	for( int nm=0; nm<_PSQH_; nm++ )
	{
		for( int ls=0; ls<_PSQH_; ls++ )
		{
			out = Y[nm*_PSQH_+ls];
			// out = Y[ IDK( l,s ) + id[n] + m];
			if( fabs(out )<1e-5) out = 0.0;
			
			cout <<out<<" ";
		} //ls
		cout<<endl;
	}// nm
}	// end printYFull

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::testY( double *Y1, double *Y2 ) const
{
	double out;
	for( int n=0; n<N_POLES; n++ )
	{
		for( int m=0; m<=n; m++ )
		{
			for( int l=0; l<N_POLES; l++ )
			{
				for( int s=0; s<=l; s++ )
				{
					if( n > l || (n==l  && m > s )) out = 0.0;
					else
					{
						out = Y1[ IDK( l,s ) + id[n] + m] + Y2[ IDK(l,s) + id[n] + m];;
						if( fabs(out )<1e-5) out = 0.0;
					}
					cout <<out<<" ";
				} // s
			} // l
			cout<<endl;
		}// m
	}// n
}	// end testY

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::testYFull( double *Y1, double *Y2 ) const
{
	int nm = 0;
	for( int n=0; n< N_POLES; n++ )
	{
		for( int m=0; m<=n; m++, nm++ )
		{
			int k = nm*_PSQH_;
			cout <<n<<","<<m<<": ";
			for( int l=0; l<N_POLES; l++ )
			{
				for( int s=0; s<=l; s++, k++ )
				{
					double out = Y1[k] + Y2[k];
					if(  !(l==n && m==s ) )
						if (  fabs(out )>1e-4)
							cout <<"( "<<l<<","<<s<<" ) "<<out<<" ";
				} // s
			} // l
			cout<<endl;
		}// m
	}
}	// end testYFull

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::printMat( double *mat )
{
	for( int i=0; i<_PSQ_; i++ )
	{
		for( int j=0; j<_PSQ_; j++ )
		{
			int k = j*_PSQ_+i;
			double out = mat[k];
			if( fabs(out ) < 1e-5) out = 0.0;
			cout <<out<<" ";
		}
		cout<<endl;
	}
	return;
}	// end printMat

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::printMat( const vector<double> &mat )
{
	int len = ( int ) sqrt( double(mat.size()));
	
	for( int i=0; i<len; i++ )
	{
		for( int j=0; j<len; j++ )
		{
			int k = j*len+i;
			double out = mat[k];
			if( fabs(out ) < 1e-12) out = 0.0;
			cout <<out<<" ";
		}
		cout<<endl;
	}
	
	return;
}	// end printMat

/*#########################################################*/
/*#########################################################*/
// debugging
/*#########################################################*/
/*#########################################################*/

#ifdef __DEBUG__

/*#########################################################*/
/*#########################################################*/
// DEBUG: : check results using Dirichlet BC
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::checkSolveDifferenceDirichlet(  )
{
	
	assert( m_kappa ==0 ); // for now the FACT constants are hardcoded
	const char* fname;
	
	if( !m_flag )
	{
		if( m_id == 0 )
			fname = "diffD_sp0_cy1.dat";
		else if ( m_id == 1 )
			fname = "diffD_sp1_cy1.dat";
		
		ofstream dfile( fname );
		
		double avgErrorE = 0.0;
		double maxErrorE = 0.0;
		double sum=0.0;
		int ct = 0;
		int npe = m_SPE.size(  );
		int npb = m_SPB.size(  );
		int freq = ( npe+npb ) / 10000;
		double delta = 0.1; // cutoff for printing
		
		for( int h=0; h<npe; h+=freq )
		{
			CSpPnt q = CartToSph( m_SPE[h] );
			CSHExpan Y( q.theta( ), q.phi(), N_POLES);
			
			double diff=0, LHS=0, RHS=0;
			for( int n=0; n<N_POLES; n++ )
			{
				LHS += Y( n,0 ) * ( m_Efix(n,0) + m_F(n,0) + m_rad
								   * (  m_Lfix(n,0 ) + m_LF(n,0) )  );
				RHS += Y( n,0 ) * ( m_H(n,0)               + m_rad * m_LH(n,0) ) ;
				
				for( int m=1; m<2*n+1; m++ )
				{
					LHS += 2 * Y( n,m ) * (m_Efix(n,m) + m_F(n,m)
										   + m_rad *(  m_Lfix(n,m ) + m_LF(n,m)) ) ;
					RHS += 2 * Y( n,m ) * (m_H(n,m)               + m_rad * m_LH(n,m) ) ;
				}
			}
			
			diff = LHS - RHS;
			sum += fabs( LHS ) + fabs(RHS);
			avgErrorE += fabs( diff );
			
			if(  fabs(diff ) > delta )
				dfile <<h<<" "<<diff<<" "<<m_SPE[h].x(  )<<" "
				<< m_SPE[h].y(  )<<" "<< m_SPE[h].z()<<endl;
			if( fabs(diff ) > maxErrorE)
				maxErrorE = fabs( diff );
			
			ct++;
		}
		
		avgErrorE /= ct;
		sum /= ct;
		
		dfile.close(  );
		
		cout <<">>>>>> SOLVE: Dirchlet AvErr "<<avgErrorE<<" MaxErr "<<maxErrorE<<
		" rel avg "<<avgErrorE/sum<<" rel max "<<maxErrorE/sum<<endl;
	}
	return;
} // end checkSolveDifferenceDirichlet

/*#########################################################*/
/*#########################################################*/
// debugging: check results using Von Neumann BC
/*#########################################################*/
/*#########################################################*/

void
CSolExpCenter::checkSolveDifferenceVonNeumann(  )
{
	const char* fname;
	
	if( !m_flag )
	{
		if( m_id == 0 )
			fname = "diffV_sp0_cy1.dat";
		else if ( m_id == 1 )
			fname = "diffV_sp1_cy1.dat";
		
		ofstream dfile( fname );
		
		double avgErrorE = 0.0;
		double maxErrorE = 0.0;
		double sum=0.0;
		int ct = 0;
		int npe = m_SPE.size(  );
		int npb = m_SPB.size(  );
		int freq = ( npe+npb ) / 10000;
		double delta = 0.1; // cutoff for printing
		
		const double eps =  m_idiel/ m_sdiel;
		double FACTE_E[N_POLES],FACTL_E[N_POLES],
		FACTH_E[N_POLES], FACTLH_E[N_POLES];
		
		assert( m_kappa ==0 );
		for( int n=0; n<N_POLES; n++ )
		{
			FACTE_E[n]  = -eps * ( n+1 );
			FACTL_E[n]  =  eps * n * m_rad;
			FACTH_E[n]  = -( n+1 );
			FACTLH_E[n] = m_rad * n;
		}
		
		for( int h=0; h<npe; h+=freq )
		{
			CSpPnt q = CartToSph( m_SPE[h] );
			CSHExpan Y( q.theta( ), q.phi(), N_POLES);
			
			double diff=0, LHS=0, RHS=0;
			for( int n=0; n<N_POLES; n++ )
			{
				LHS += Y( n,0 ) * ( FACTE_E[n]*m_Efix(n,0) + eps*n*m_F(n,0) +
								   FACTL_E[n]*(  m_Lfix(n,0 ) + m_LF(n,0))  );
				RHS += Y( n,0 ) * ( FACTH_E[n]*m_H(n,0) + FACTLH_E[n]*m_LH(n,0) ) ;
				
				for( int m=1; m<2*n+1; m++ )
				{
					LHS += 2 * Y( n,m ) * (FACTE_E[n]*m_Efix(n,m) + eps*n*m_F(n,m) +
										   FACTL_E[n]*(  m_Lfix(n,m ) + m_LF(n,m))) ;
					RHS += 2 * Y( n,m ) * (FACTH_E[n]*m_H(n,m) + FACTLH_E[n]*m_LH(n,m) ) ;
				}
			}
			diff = LHS - RHS;
			sum += fabs( LHS ) + fabs(RHS);
			avgErrorE += fabs( diff );
			
			if(  fabs(diff ) > delta )
				dfile <<h<<" "<<diff<<" "<<m_SPE[h].x(  )<<" "
				<< m_SPE[h].y(  )<<" "<< m_SPE[h].z()<<endl;
			if( fabs(diff ) > maxErrorE)
				maxErrorE = fabs( diff );
			ct++;
		} // end h
		avgErrorE /= ct;
		sum /= ct;
		
		dfile.close(  );
		cout <<">>>>>> SOLVE: Von Neumann AvErrE "<<avgErrorE<<" MaxErrE "<<maxErrorE<< 
		" rel avg "<<avgErrorE/sum<<" rel max "<<maxErrorE/sum<<endl;
	}

	return;
}	// end checkSolveDifferenceVonNeumann

#endif // __DEBUG

