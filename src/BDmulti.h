#ifndef _BDM_H_
#define _BDM_H_
#include <vector>
#include "util.h"
#include "readutil.h"
#include "molecule.h"
#include "system.h"

#define FACT2 5.604586e2 //( COUL_K*IKbT )
const double AVOGADRO_NUM = 6.0221415e23;
const double ELECT_CHG = 1.60217646e-19;
const double PERMITTIVITY_VAC = 8.854187817e-12;
const double KB = 1.3806503e-23;
const double ANGSTROM = 1e-10;
const double LITRE = 1e-3;
const double IKbT = 1.0/( KB*298.15 );
const double COUL_K = 2.30708e-18; // ( 1/4PI*e0 ) * e2 / ANGSTROM

using namespace std;

//!  The CBDmulti class
/*! The class that contains information about a BD multi run */
class CBDmulti
{
public:
	// Initialization
	CBDmulti( int Np1, int Np2, const vector<char*> &molfnames1,
					 const vector<char*> &molfnames2, REAL idiel );
	
	static void initConstants( int nMolType );
	static void addMolContact( int mol1type,
														const vector<CMolContact>&molcontactlist );
	
	void resetLattice(  );
	void restartConfig( char* configname );
	
	double runFirstDockTime( double maxtime, char* fname );
	CMolecule* getMol( int i ) {return m_mols[i];}
	static CPnt getRandVec( REAL std )
	{ return std*CPnt( normRand( ), normRand(), normRand()); }
	
	vector<CMolecule*> m_mols; // debug in public
	vector<vector<CPnt> > m_SPxes;
	vector<CMolCell> m_molcells;
	
private:
	bool IsDocked( CMolecule* mol1, CMolecule* mol2 ) const;
	
	void initMolTypes( const vector<char*> &molfnames1,
										const vector<char*> &molfnames2, REAL idiel );
	
	double propagate( double mindist, bool bWrite ); // debug
	double propagate( double mindist ) { return propagate(mindist, false);}
	
	REAL computeMinSepDist(  );
	REAL compute_dt( double minsepdist ) const;
	void makeMove( const vector<CPnt> & dR, const vector<CPnt> & dO, REAL dt );
	
	static REAL MINDIST, MAXDIST;
	
	static int m_nMolType;
	static vector< vector<CMolContact> > MOLCONTACTLIST;
	
	// moltype variables
	vector<REAL*> m_iMats1;												//!< A matrix of interaction matrices for molecule1
	vector<REAL*> m_iMats2;												//!< A matrix of interaction matrices for molecule2
	vector<CMulExpan*> m_iF1;											//!< A multipole expansion for F of molecule1
	vector<CMulExpan*> m_iH1;											//!< A multipole expansion for H of molecule1
	vector<CMulExpan*> m_iF2;											//!< A multipole expansion for F of molecule2
	vector<CMulExpan*> m_iH2;											//!< A multipole expansion for H of molecule2
	vector<vector<double> > m_qSolvedF1;					//!< A vector of solved F for molecule1
	vector<vector<double> > m_qSolvedH1;					//!< A vector of solved H for molecule1
	vector<vector<double> > m_qSolvedF2;					//!< A vector of solved F for molecule2
	vector<vector<double> > m_qSolvedH2;					//!< A vector of solved H for molecule2
	vector<double> m_totalF1;											//!< A vector of summed F for molecule1
	vector<double> m_totalH1;											//!< A vector of summed H for molecule1
	vector<double> m_totalF2;											//!< A vector of summed F for molecule2
	vector<double> m_totalH2;											//!< A vector of summed H for molecule2
	vector<CLocalExpan*> m_LFs_intraSelf1;				//!< A local expansion for F of molecule1
	vector<CLocalExpan*> m_LHs_intraSelf1;				//!< A local expansion for H of molecule1
	vector<CLocalExpan*> m_LFs_intraSelf_far1;		//!< A local expansion for far F of molecule1
	vector<CLocalExpan*> m_LHs_intraSelf_far1;		//!< A local expansion for far H of molecule1
	vector<CLocalExpan*> m_LFs_intraSelf2;				//!< A local expansion for F of molecule2
	vector<CLocalExpan*> m_LHs_intraSelf2;				//!< A local expansion for H of molecule2
	vector<CLocalExpan*> m_LFs_intraSelf_far2;		//!< A local expansion for far F of molecule2
	vector<CLocalExpan*> m_LHs_intraSelf_far2;		//!< A local expansion for far H of molecule2
	vector<CMolCell> m_molcells1;									//!< A vector of cells for molecule1
	vector<CMolCell> m_molcells2;									//!< A vector of cells for molecule2
	vector<vector<int> > m_neighs1;								//!< A vector of neighbors for each sphere in molecule1
	vector<vector<int> > m_neighs2;								//!< A vector of neighbors for each sphere in molecule2
	vector<vector<CPnt> > m_SPxes1;								//!< A vector of coordinates of exposed surface points in molecule1
	vector<vector<CPnt> > m_SPxes2;								//!< A vector of coordinates of exposed surface points in molecule2
	vector<int> m_nSPx1;													//!< A vector of indices of exposed surface points in molecule1
	vector<int> m_nSPx2;													//!< A vector of indices of exposed surface points in molecule2
	vector< vector<int> > m_intraPolLists_near1;	//!< A vector of lists of spheres considered near in molecule1
	vector< vector<int> > m_intraPolLists_near2;	//!< A vector of lists of spheres considered near in molecule2
	
	// BD variables
	int m_np1;						//!< The number of molecules of type 1
	int m_np2;						//!< The number of molecules of type 2
	vector<REAL> m_Dtr;		//!< The translational diffusion coefficient
	vector<REAL> m_Dr;		//!< The rotational diffusion coefficient
	vector<CQuat> m_rots;	//!< The rotational orientation
	CPnt m_initialrcen1;	//!< The xyz coord of initial read-in of molecule 1
	CPnt m_initialrcen2;	//!< The xyz coord of initial read-in of molecule 1
}; // end CBDmulti

#endif
