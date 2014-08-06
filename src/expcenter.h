#ifndef _EXPCENTER_H_
#define _EXPCENTER_H_
#include "expansion.h"
#include "rotcoeff.h"
#include "hash.h"

#define _PSQH_ N_POLES*( N_POLES+1 )/2			//!< The result of (NPOL*(NPOL+1))/2
#define _PSQH2_ _PSQH_ * _PSQH_					//!< The result of ((NPOL*(NPOL+1))/2)^2
#define _PQUADH_ _PSQH_ * ( _PSQH_+1 ) / 2		//!< The result of (NPOL*(NPOL+1)^2)/4
#define _PSQ_ N_POLES*N_POLES						//!< The result of NPOL^2
#define _PQUAD_ _PSQ_*_PSQ_							//!< The result of NPOL^4

//! CExpCenter class
/*!	This class is the base class for all expansion centers */
class CExpCenter
{
public:
	//! CExpCenter initConstants function
	/*!	This function is used to initialize kappa and sdiel. 
	 Also calls an initialization of the RExpan object constants
	 \param sdiel a floating point of the solvent's dielectric constant
	 \param kappa a floating point of the inverse debye length */
	static void initConstants( double sdiel, double kappa );
	
	CExpCenter(  );
	CExpCenter( int ki, CPnt cen, double rad );
	
	CExpCenter( int ki, CPnt cen, double rad, double idiel, bool bKappa,
			   const vector<double> &chgAssigned,
			   const vector<CPnt> &posAssigned,
			   CQuat* orient=NULL, CRotCoeff* rot=NULL );
	
	void setOrder( const int p, const CRotCoeff &rot );
	void incOrder( const CRotCoeff & rot );
	const REAL computePot(  ) const;
	
	const double getRad(  ) const { return m_rad;}
	const CPnt getCen(  ) const { return m_cen;}
	const CPnt getCenRot(  ) const { return m_cenRot;}
	const CQuat getOrient(  ) const {return *m_orient;}
	const int getOrder(  ) const {return m_p;}
	const double getLScale(  ) const {return m_lscale;}
	const double getMScale(  ) const {return m_mscale;}
	const bool getbKappa(  ) const {return m_bKappa;}
	
	const CMulExpan * const getpH(  ) const  {return &m_H;}
	CMulExpan * const getpH(  ) {return &m_H;}
	CMulExpan & getH(  ) {return m_H;}
	const CMulExpan getH(  ) const {return m_H;}
	CMulExpan & getHrot(  ) {return m_Hrot;}
	const CMulExpan getHrot(  ) const {return m_Hrot;}
	const double getHrotMono(  ) const {return m_Hrot[0];}

	CLocalExpan & getLS(  ) {return m_LS;}
	const CLocalExpan & getLS(  ) const {return m_LS;}
	const CMulExpan getrE(  ) const {return m_rE;}
	
	CGradExpan & getgLHN(  ) {return m_gLHN;}
	const CGradExpan getgLHN(  ) const {return m_gLHN;}
	
	CExpCenter & operator=( const CExpCenter &M );
	
	// debug
	const CRotCoeff getRotCoeff(  ) const {return *m_rot;}
	
protected:
	void setFixedCharges(  );
	
	static double m_sdiel;			//!< A floating point of the solvent dielectric constant
	static double m_kappa;			//!< A floating point of the system's inverse debye length
	CPnt m_cen;						//!< An xyz of the center of the molecule in the  standard frame
	CPnt m_cenRot;					//!<
	double m_rad;						//!< A floating point of the current CG sphere's radius (a^(I,k))
	double m_lscale;					//!< A floating point of the length scale of the current sphere, defaulted to m_rad
	double m_mscale;					//!< A floating point of another length scale, defaulted to 1/m_rad
	double m_idiel;					//!< A floating point of the molecule's interior dielectric constant
	bool m_bKappa;					//!<
	vector<double> m_ch;				//!<
	vector<CPnt> m_pos;				//!<
	
	CLocalExpan m_LS;				//!< only due to external source ( selected )
	CMulExpan m_E;					//!<
	CMulExpan m_rE;					//!<
	CMulExpan m_H;					//!< A multipole expansion of the effective polarization charges on sphere k of the exposed surf
	CMulExpan m_Hrot;				//!< A rotated MPE of the effective polarization charges on sphere k of the exposed surf
	CGradExpan m_gH;					//!<
	CGradExpan m_gLHN;				//!<
	
	CQuat* m_orient;					//!< A quaterion of the the current rotation of the molecule
	CRotCoeff* m_rot;				//!< 
	
	int m_p;							//!< An integer of the number of poles in the expansion
	int m_id;							//!< An integer of the ID number of the CG sphere
};  // end  class CExpCenter


//! CSolExpCenter class
/*!	This class is for expansion centers for molecules with partially exposed
 surfaces  */
class CSolExpCenter : public CExpCenter
{
public:
	//! CSolExpCenter class constructor
	/*!	This constructor is used for generating IMATRICES
	 \param ki an integer of the CG sphere of interest
	 \param cen an XYZ vector of the sphere of interest
	 \param rad a floating point of radius the spbere of interest
	 \param SPx a vector of the expansion points for sphere ki
	 \param nSPx an integer of the number of expansion points for ki
	 \param neigh a vector of integers indicating the IDs of the neighboring CG spheres
	 \param intraPolList_near a list of the CG spheres considered near
	 \param idiel a floating point of the interior dielectric
	 \param bKappa a boolean of whether or not an inverse debye length is given
	 \param allchg a vector of all floating point charges in the molecule
	 \param allpos a vector of XYZ positions of point charges WRT to ki CG sphere
	 \param chgAssigned	a vector of all floating point charges within sphere ki
	 \param posAssigned	a vector of coords of all floating point charges within sphere ki  */
	CSolExpCenter( int ki, CPnt cen, double rad, const vector<CPnt> &SPE,
				  const vector<CPnt> &SPB,
				  const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
				  const vector<int> &intraPolList_near,
				  double idiel, bool bKappa, const vector<double> &allchg,
				  const vector<CPnt> &allpos,
				  const vector<double> &chgAssigned, const vector<CPnt> &posAssigned );
	//! CSolExpCenter class constructor
	/*!	This constructor is used for generating SELFPOL
	 \param ki an integer of the CG sphere of interest
	 \param cen an XYZ vector of the sphere of interest
	 \param rad a floating point of radius the spbere of interest
	 \param SPx a vector of the expansion points for sphere ki
	 \param nSPx an integer of the number of expansion points for ki
	 \param neigh a vector of integers indicating the IDs of the neighboring CG spheres
	 \param intraPolList_near a list of the CG spheres considered near
	 \param idiel a floating point of the interior dielectric
	 \param bKappa a boolean of whether or not an inverse debye length is given
	 \param allchg a vector of all floating point charges in the molecule
	 \param allpos a vector of XYZ positions of point charges WRT to ki CG sphere
	 \param chgAssigned	a vector of all floating point charges within sphere ki
	 \param posAssigned	a vector of coords of all floating point charges within sphere ki 
	 \param iMat a vector of Imatrices for setting up the self polarization  */
	CSolExpCenter( int ki, CPnt cen, double rad,
				  const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
				  const vector<int> &intraPolList_near,
				  double idiel, bool bKappa, const vector<double> &allchg,
				  const vector<CPnt> &allpos,
				  const vector<double> &chgAssigned, const vector<CPnt> &posAssigned,
				  REAL* iMat );	
	//! CSolExpCenter class constructor
	/*!	This constructor
	 \param ki an integer of the CG sphere of interest
	 \param cen an XYZ vector of the sphere of interest
	 \param rad a floating point of radius the spbere of interest
	 \param SPx a vector of the expansion points for sphere ki
	 \param nSPx an integer of the number of expansion points for ki
	 \param neigh a vector of integers indicating the IDs of the neighboring CG spheres
	 \param intraPolList_near a list of the CG spheres considered near
	 \param idiel a floating point of the interior dielectric
	 \param bKappa a boolean of whether or not an inverse debye length is given
	 \param allchg a vector of all floating point charges in the molecule
	 \param allpos a vector of XYZ positions of point charges WRT to ki CG sphere
	 \param chgAssigned	a vector of all floating point charges within sphere ki
	 \param posAssigned	a vector of coords of all floating point charges within sphere ki
	 \param iMat  a vector of Imatrices
	 \param Fself an object of coefficients solved for with self polarization
	 \param Hself an object of coefficients solved for with self polarization
	 \param LF_intraSelf a local expansion object from self pol
	 \param LH_intraSelf a local expansion object from self pol
	 \param LF_intraSelf_far a local expansion object from self pol
	 \param LH_intraSelf_far a local expansion object from self pol
	 \param qSolvedFself a floating point of q from F
	 \param qSolvedHself a floating point of q from H
	 \param totalFself a floating point of the sum of F self-pol coefficients
	 \param totalHself a floating point of the sum of H self-pol coefficients
	 \param orient a quaternion of the expansion's orientation
	 \param rot a RotCoeff object of the expansion's rotation coefficients
	 \param bGrad a boolean indicating whether or not to compute the gradient */
	CSolExpCenter( int ki, CPnt cen, double rad,
				  const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
				  const vector<int> &intraPolList_near,
				  double idiel, bool bKappa, const vector<double> &allchg,
				  const vector<CPnt> &allpos,
				  const vector<double> &chgAssigned,
				  const vector<CPnt> &posAssigned, REAL* iMat,
				  const CMulExpan *Fself, const CMulExpan *Hself,
				  const CLocalExpan* LF_intraSelf, const CLocalExpan* LH_intraSelf,
				  const CLocalExpan* LF_intraSelf_far,
				  const CLocalExpan* LH_intraSelf_far,
				  const double *qSolvedFself, const double  *qSolvedHself,
				  const double *totalFself, const double  *totalHself,
				  CQuat* orient, CRotCoeff* rot, bool bGrad );
	
	CSolExpCenter( int ki, CPnt cen, double rad,
				  const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
				  const vector<int> &intraPolList_near,
				  const CMulExpan *Fself, const CMulExpan *Hself,
				  const double *qSolvedFself, const double  *qSolvedHself,
				  const double *totalFself, const double  *totalHself );
	
	// for queries
	CSolExpCenter( int ki, CPnt cen, double rad, const CMulExpan *Hself,
				  const vector<CPnt> &SPx=vector<CPnt>( 0 ),
				  const int &nSPx=0, const vector<int> & intraPolList_near=vector<int>( 0 ));
	
	~CSolExpCenter(  )
	{ if( m_bOwnMat ) delete [] m_IMat; }
	
	// public functions
	
	//! CSolExpCenter initConstants function
	/*!	This function is used to initialize many system constants and conditions,
	 like (2n+1)/(4*pi) in EQ 15-17 in Yap 2010.
	 \param sdiel a floating point of the solvent's dielectric constant
	 \param kappa a floating point of the inverse debye length */
	static void initConstants( double sdiel, double kappa );
	static void build_IMat( const vector<CPnt> &SPE, const vector<CPnt> &SPB, REAL* IMat );
	static void computeExposedSurfaceChargesFH_( const vector<CPnt> &SPx,
												const CMulExpan & F, const CMulExpan & H,
												double *constH, double nSPx,
												vector<double> &qSolvedF, vector<double> &qSolvedH,
												double &totalF, double &totalH, int p );
	
	//! CSolExpCenter solveSurfaceCharges function
	/*!	This function is */
	double solveSurfaceCharges(  );
	double solveSurfaceCharges( const CMulExpan &Finit,
							   const CMulExpan &Hinit, bool bRotateH );
	double solveSurfaceGradient( int j, CGradExpan &gF, CGradExpan &gH,
								CGradExpan &gHrot,
								const CGradExpan & LGF_j, const CGradExpan & LGH_j );
	void initGradient( int nmol );
	void resetExpansions( const CMulExpan &iF, const CMulExpan &iH );
	void setToSelfPolValues( bool bResetCharges );
	const vector<REAL> & getExposedSurfaceChargesH(  ) const
	{return m_qSolvedH; }
	const vector<REAL> & getExposedSurfaceChargesF(  ) const
	{return m_qSolvedF; }
	
	void rotateCenters(  );
	void rotateHself(  );
	
#ifdef __NOPOL_UNCHANGED__
	void rotateCurrentH(  );
	CGradExpan rotateGH(  ) const;
#endif
	
	REAL getDev(  ) const {return m_dev;}
	CMulExpan & getF(  ) {return m_F;}
	const CMulExpan getF(  ) const {return m_F;}
	
	CGradExpan & getGH(  ) {return m_gH;}
	const CGradExpan getGH(  ) const {return m_gH;}
	
	const CPnt getGHMono(  ) const
	{return CPnt( m_gH[0][0],m_gH[1][0],m_gH[2][0] );}
	
	CLocalExpan & getLHS(  ) {return m_LHS;}
	const CLocalExpan getLHS(  ) const {return m_LHS;}
	CLocalExpan & getLFS(  ) {return m_LFS;}
	const CLocalExpan getLFS(  ) const {return m_LFS;}
	
	CLocalExpan & getLFS_Far(  ) {return m_LFS_Far;}
	const CLocalExpan getLFS_Far(  ) const {return m_LFS_Far;}
	CLocalExpan & getLHS_Far(  ) {return m_LHS_Far;}
	const CLocalExpan getLHS_Far(  ) const {return m_LHS_Far;}
	
	const vector<CPnt> & getSPx(  ) const {return m_SPx;}
	
	const vector<double> & getQSolvedF(  ) const {return m_qSolvedF;}
	const vector<double> & getQSolvedH(  ) const {return m_qSolvedH;}
	const vector<CPnt> & getGSolvedF(  ) const {return m_gSolvedF;}
	const vector<CPnt> & getGSolvedH(  ) const {return m_gSolvedH;}
	const int getSPExSize(  ) const {return m_SPExSize;}
	
	const REAL& getTQH(  ) const {return m_totalH;}
	const REAL& getTQF(  ) const {return m_totalF;}
	const REAL& getTGH(  ) const {return m_maxgH;}
	const REAL& getTGF(  ) const {return m_maxgF;}
	
	const CMulExpan &getFself(  ) const {return *m_Fself;}
	const CMulExpan &getHself(  ) const {return *m_Hself;}
	const CMulExpan &getHselfRot(  ) const {return m_HselfRot;}
	const CLocalExpan getLF_intraSelf(  ) const {return *m_LF_intraSelf;}
	const CLocalExpan getLH_intraSelf(  ) const {return *m_LH_intraSelf;}
	const CLocalExpan getLF_intraSelf_far(  ) const {return *m_LF_intraSelf_far;}
	const CLocalExpan getLH_intraSelf_far(  ) const {return *m_LH_intraSelf_far;}
	
	const REAL* getQSolvedFself(  ) const {return m_qSolvedFself;}
	const REAL* getQSolvedHself(  ) const {return m_qSolvedHself;}
	const REAL getTQFself(  ) const {return *m_totalFself;}
	const REAL getTQHself(  ) const {return *m_totalHself;}
	
	CSolExpCenter & operator=( const CSolExpCenter &M );
	
	const vector<int> & getIntraPolList_near(  )
	{return m_intraPolList_near;}
	
	void clearInterPolList(  )
	{m_interactList.clear(  ); m_interPolList.clear();	}
	void addInterPolList( CFullSphereID c ) {m_interPolList.push_back(c);}
	void addInteractList( CFullSphereID c ) {m_interactList.push_back(c);}

	const vector<CFullSphereID> & getInterPolList(  ) const
	{return m_interPolList;}
	const vector<CFullSphereID> & getInteractList(  ) const
	{return m_interactList;}

	
#ifdef __NOPOL_UNCHANGED__
	void setbInterPolListChanged( bool bChanged )
	{m_bInterPolListChanged = bChanged;}
	bool IsInterPolListChanged(  ) const
	{return m_bInterPolListChanged;}
#endif
	
	bool IsEmptyInterPolList(  ) const { return m_interPolList.size() == 0;}
	bool IsNoInteractionList(  ) const { return m_interPolList.size() == 0 &&
		m_interactList.size(  ) == 0;  }
	
	bool isOnInterPolList( int j, int kj ) const;
	
	virtual const CPnt computeTorqueOn_0( int i ) const;
	const CPnt computeForceOn_0(  ) const;
	
	// public variables!
	static double CONST2[N_POLES];						//!< The constant (2l+1)/(4*pi) used in EQs 15-17 in Yap 2010
	
private:
	// private functions
	static int & IDK( int l, int s )    {return idk[id[l]+s];}

	//! CSolExpCenter computeIntegralE function
	/*!	This function is used to compute surface integra; directly using
	 the exposed surface points.
	 \param Yrr a vector of the real components of the bessel multiplication
	 \param Yri a vector of the real of bessel ls * imag(Conj nm)
	 \param Yir a vector of the imag of bessel ls * real(Conj nm)
 	 \param Yii a vector of the imag of bessel ls * imag(Conj nm)
	 \param SPE the exposed surface points 
 	 \param SPB the buried surface points */
	static void computeIntegralE( double * Yrr,double * Yri,
								 double * Yir,double * Yii,
								 const vector<CPnt> & SPE, const vector<CPnt> & SPB );
	//! CSolExpCenter computeIntegralB function
	/*!	This function is used to compute surface integra; indirectly using
	 the buried surface points.
	 \param Yrr a vector of the real components of the bessel multiplication
	 \param Yri a vector of the real of bessel ls * imag(Conj nm)
	 \param Yir a vector of the imag of bessel ls * real(Conj nm)
 	 \param Yii a vector of the imag of bessel ls * imag(Conj nm)
	 \param SPE the exposed surface points
 	 \param SPB the buried surface points */
	static void computeIntegralB( double * Yrr,double * Yri,
								 double * Yir,double * Yii,
								 const vector<CPnt> & SPE, const vector<CPnt> & SPB  );
	//! CSolExpCenter setFixedCharges function
	/*!	This function is used to initialize many system constants for EQ 22 in Yap 2010.  
	 It calls CMul and Local expansions to create them for the system. 
	 \param allchg a vector of all the charges in the molecule
	 \param a vector of coordinates for each point charge of the molecule, in ref
				frame of current	CG sphere */
	virtual void setFixedCharges( const vector<double> &allchg,
								 const vector<CPnt> &allpos );
	
	//! CSolExpCenter initMyConstants function
	/*!	This function is used to for a given CG sphere, 1. generate m_E, m_H
	 2. set charges for current center and generate m_Efix, m_Lfix.   */
	void initMyConstants(  );
	//! CSolExpCenter initSurfaceCharges function
	/*!	This function is used to initialize surface charges for polarization   */
	void initSurfaceCharges(  );
	//! CSolExpCenter computeExposedSurfaceChargesH function
	/*!	This function is used to initialize H surface charges for polarization   */
	void computeExposedSurfaceChargesH(  );
	//! CSolExpCenter computeExposedSurfaceChargesF function
	/*!	This function is used to initialize F surface charges for polarization   */
	void computeExposedSurfaceChargesF(  );
	void computeExposedSurfaceChargesFH(  );
	void computeExposedSurfaceGradientFH( const CGradExpan &gF,
										 const CGradExpan &gH );
	//! CSolExpCenter build_IMatE function
	/*!	This function is used to build a response matrix for exposed surface points  */
	void build_IMatE(  );
	
	void cal_DVFixed(  );
	void getSurfaceChargesH( const vector<CPnt> &SPx, vector<double> &qS );
	void getSurfaceChargesF( const vector<CPnt> &SPx, vector<double> &qS );

	//! CSolExpCenter revertF function
	/*!	This function reverts the F matrix by multiplying it by (1/(2l+1))
	 \param F a vector of F constants
	 \param p the number of poles
	 \param D an int for the number of dimensions */
	void revertF( double * F, int p, int D ) const ;
	//! CSolExpCenter revertF function
	/*!	This function reverts the H matrix by multiplying it by (ihat_n(kr)/(2n+1))
	 \param H a vector of H constants
	 \param p the number of poles
	 \param D an int for the number of dimensions */
	void revertH( double * H, int p, int D ) const ;
	
	double solveFH_N( const vector<double> &LH,
					 const vector<double>  &LF, int pm, int pn );
	//! CSolExpCenter solveFH function
	/*!	This function iteratively solves for F, H using D, V separately
	 conjugation taken into account in IMAT
	 I think this is the core function for doing mutual polarization
	 In general, what function 'solveFH' did is to calculate F and H
	 iteratively with given input of F, H, LF and LH.
	 It has called function applyMM, which can be found in file 'lautil.cpp',
	 to do matrix vector products. So the job 'solveFH' has done is to
	 solve linear equations in 4a and 4b in 2013 JCTC paper of THG and Enghui Yap.
	 
	 Note that 'solveFH' is not only used for solving multiple expansion of
	 surface charge density like F and H, but also the gradients of multiple
	 Add a comment to this line expansions, such as equation 11a and 11b in 2013
	 JCTC paper of THG and Enghui Yap. ( S. Liu )
	 \param F a vector of F constants
	 \param H a vector of H constants
	 \param LF a vector of local f constants
	 \param LH a vector of local H constants
	 \param pm the number of poles 
	 \param pn the number of poles
	 \param D an int for the number of dimensions */
	double solveFH( vector<double> &F, vector<double>  &H,
				   const vector<double> &LF, const vector<double>  &LH,
				   int pm, int pn, int D ) const;
	//! CSolExpCenter compute_FHbase function
	/*!	This function computes hbase and fbase to incorporate LF, LH
	 \param LF pointer to a vector of local f constants
	 \param LH pointer to a vector of local H constants
	 \param fbase pointer to a vector of base local f constants
	 \param hbase pointer to a vector of base local H constants
	 \param pm the number of poles  */
	void compute_FHbase( const double* LH, const double *LF, double* fbase,
						double *hbase, int pm ) const;
	//! CSolExpCenter compute_FHbase3 function
	/*!	This function computes hbase and fbase to incorporate LF, LH
	 \param LF pointer to a vector of local f constants
	 \param LH pointer to a vector of local H constants
	 \param fbase pointer to a vector of base local f constants
	 \param hbase pointer to a vector of base local H constants
	 \param pm the number of poles  */
	void compute_FHbase3( const double* LH, const double *LF, double* fbase,
						 double *hbase, int pm ) const;
	//! CSolExpCenter compute_fx function
	/*!	This function computes fbase to incorporate into fx
	 \param fx pointer to a vector of f expansion constants
	 \param fbase pointer to a vector of base local f constants
	 \param F pointer to a vector of f constants	
	 \param H pointer to a vector of h constants	 
	 \param hbase pointer to a vector of base local H constants
	 \param p the number of poles 
	 \param D an int for the number of dimensions  */
	void compute_fx( double *fx, const double *fbase, const double *F,
					const double *H, int p, int D ) const;
	//! CSolExpCenter compute_hx function
	/*!	This function computes hbase to incorporate into hx
	 \param hx pointer to a vector of h expansion constants
	 \param hbase pointer to a vector of base local h constants
	 \param F pointer to a vector of f constants
	 \param H pointer to a vector of h constants
	 \param hbase pointer to a vector of base local H constants
	 \param p the number of poles
	 \param D an int for the number of dimensions  */
	void compute_hx( double *hx, const double *hbase, const double *F,
					const double *H, int p, int D ) const;
	//! CSolExpCenter computeDev function
	/*!	This function computes the deviation between 2 matrices
	 \param M1 a pointer to matrix 1
	 \param M2 a pointer to matrix 2
	 \param p the number of poles
	 \param D an int for the number of dimensions  */
	static double computeDev( const double *M1, const double *M2, int p, int D );
	//! CSolExpCenter computeDev2 function
	/*!	This function computes the deviation between 2 matrices
	 \param M1 a pointer to matrix 1
	 \param M2 a pointer to matrix 2
	 \param p the number of poles
	 \param D an int for the number of dimensions  */
	static double computeDev2( const double *M1, const double *M2, int p, int D );
	
	// utilities
	static void printY( double *Y ) ;
	static void printYFull( double *Y ) ;
	static void printMat( double *mat ) ;
	static void printMat( const vector<double> &mat ) ;
	void testY( double *Y1, double *Y2 ) const;
	void testYFull( double *Y1, double *Y2 ) const;
	void checkSolveDifferenceDirichlet(  );
	void checkSolveDifferenceVonNeumann(  );
	
	static double & getYY( double * YY, int n, int m, int l, int s )
	{return YY[ IDK( l,s ) + id[n] + m];}
	
	// private variables
	static double CONST3[N_POLES];				//!< 1/(2l+1)
	static int id[N_POLES];						//!< l*(l+1)/2
	static int idk[_PSQH_];						//!< stores index of the first row of every ( l,s ) column
	
	vector<CPnt> m_SPE;							//!< A vector of positions of surface points that are solvent exposed
	vector<CPnt> m_SPB;							//!< A vector of positions of sphere surface points that are buried
	vector<int> m_neigh;
	const vector<int> & m_intraPolList_near;
	vector<CFullSphereID> m_interPolList;
	vector<CFullSphereID> m_interactList; 
	vector<double> m_qSolvedF;
	vector<double> m_qSolvedH;
	vector<CPnt> m_gSolvedF;
	vector<CPnt> m_gSolvedH;
	REAL* m_IMat;									//!< The integral matrix computed using the quadrature method in IMat
	double m_Dfix[_PSQ_];						//!< A vector of coefficients for fixed constants in EQ 14a
	double m_Vfix[_PSQ_];						//!< A vector of coefficients for fixed constants in EQ 14b
	REAL m_totalH;								//!< The sum of the H values 
	REAL m_totalF;
	REAL m_maxgH;
	REAL m_maxgF;
	REAL m_dev;
	REAL m_dA;										//!< A ratio of (4*pi)/nSPX (amount of surface area per surface charge)
	
	double m_constH1[N_POLES];					//!< H-coeff for F-eqn in 22b Yap 2010
	double m_constF1[N_POLES];					//!< F-coeff for F-eqn in 22b Yap 2010
	double m_constH2[N_POLES];					//!< H-coeff for H-eqn in 22a Yap 2010
	double m_constF2[N_POLES];					//!< F-coeff for H-eqn in 22a Yap 2010
	double m_constLH1[N_POLES];					//!< for F-eqn, LH of XF component in 14b Yap 2010
	double m_constLF1[N_POLES];					//!< for F-eqn, LF of XF component in 14b Yap 2010
	double m_constLH2[N_POLES];					//!< for H-eqn, LH of XH component in 14a Yap 2010
	double m_constLF2[N_POLES];					//!< for H-eqn, LF of XH component in 14a Yap 2010
	double m_constInvH[N_POLES];				//!< Inverse of common constant: (ihat_n(kr)/(2n+1))
	double m_constChgH[N_POLES];				//!< Common constant in many H equs: (2n+1)/(4pi*ihat_n(kr))
	
	CLocalExpan m_Lfix;							//!< A Local expansion for point charges outside the CG sphere
	CLocalExpan m_LHS;							//!<
	CLocalExpan m_LFS;							//!< A Local expansion of reactive charges on the surface of sphere ki
	CLocalExpan m_LFS_Far;						//!<
	CLocalExpan m_LHS_Far;						//!<
	CMulExpan m_Efix;							//!< A multipole expansion for point charges outside the CG sphere
	CMulExpan m_F;								//!< A multipole expansion
	CMulExpan m_Ev;								//!<
	CMulExpan m_HselfRot;						//!<
	bool m_bRead;									//!< A boolean of whether or not to read in exposed surface charges for F
	bool m_bOwnMat;								//!<
	bool m_bGrad;									//!<
	
	const CMulExpan *m_Fself, *m_Hself;
	const CLocalExpan *m_LF_intraSelf, *m_LH_intraSelf,
	*m_LF_intraSelf_far, *m_LH_intraSelf_far;
	const double *m_qSolvedFself, *m_qSolvedHself;
	const double *m_totalFself, *m_totalHself;
	const vector<CPnt> &m_SPx;
	const int m_SPExSize;					//!< The size of the surface point expansion vector for sphere ki
	const int &m_nSPx;						//!< The size of the surface point expansion vector for sphere ki			
#ifdef __NOPOL_UNCHANGED__
	bool m_bInterPolListChanged;
#endif
};	// end class CSolExpCenter : public CExpCenter


/******************************************************************/
/******************************************************************/
/**
 * ExpCenter function setOrder a function that sets the order of the
 expansion and increases or decreases it accordingly
 ******************************************************************/
inline void
CExpCenter::setOrder( const int p, const CRotCoeff &rot )
{
	assert( p <= N_POLES );
	cout <<"setorder is used" <<endl;
	if ( m_p < p )
		while (  m_p < p )
			incOrder( rot );
	
	else if ( m_p > p )
	{
		// ( LATER ) m_rT[ki].setRange(p);
		m_H.setRange( p ); //CHECK LATER
		m_p = p;
	}
	/* ( LATER )
	 if( ki == 0 )
	 for( int k=0; k<m_ncenter; k++ )
	 m_pG[k].setRange( p );
	 */
}	// end CExpCenter::setOrder

#endif


