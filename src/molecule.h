#ifndef _MOLECULE_H_
#define _MOLECULE_H_

#include <fstream>
#include <string>
#include <sstream>
#include <cfloat>
#include <vector>
#include "expcenter.h"
#include "readutil.h"
#include "system.h"
#include "xforms.h"
#include "hash.h"
#include "contactRigid.h"
#include "contact.h"
#include "cell.h"

#ifdef __OMP
	#include "omp.h"
#endif

#define MINSEPDIST 5				//!< A minimum separation distance
#define PS_TO_NS 0.001			//!< A conversion factor: 1 ps = 0.001 ns 

class CSystem;

//!  The CMolecule class
/*!	The class that contains all details for a molecule object  */
class CMolecule
{
public:

	//! CMolecule constructor for IMatrices
	/*!	This constructor is called for IMat calculations
			\param rcen an XYZ coordinate for the center of geom of the molecule
			\param cens a vector of XYZ coordinates of CG spheres
			\param radii a vector of floating point radii of CG spheres
			\param CH a vector of floating point of point charges in the molecule
			\param POS a vector of coordinates for each point charge
			\param idiel a floating point of the dielectric of the molecule
			\param SPxes a vector of transforms for the CG spheres in the molecule
			\param nSPx a vector of indices for the transforms of the CG spheres
			\param neighs a vector of neighbor lists 
			\param intraPolLists_near a vector of spheres that are considered close
			\param molcell a vector of molecule cells - by default created empty.  */	
	CMolecule( CPnt rcen, const vector<CPnt> &cens, 
						const vector<double> &radii, const vector<double> &CH, 
						const vector<CPnt> &POS, double idiel,
						const vector<vector<CPnt> > &SPxes,  const  vector<int> &nSPx, 
						const vector<vector<int> >&neigh, 
						const vector< vector<int> > &intraPolLists_near /*, 
						const vector<CMolCell> &molcell = vector<CMolCell>( 0 )*/);
	
	CMolecule( CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii, 
						const vector<double> &CH, 
						const vector<CPnt> &POS, double idiel,  vector<REAL*> &iMats, 
						REAL intraRcutoff, 
						const vector<vector<CPnt> > &SPxes,  const  vector<int> &nSPx, 
						const vector<vector<int> >&neighs,
						const vector< vector<int> > &intraPolLists_near, 
						const vector<CMolCell> &molcell = vector<CMolCell>( 0 ));
	
	// for multi and nafion
	CMolecule( int moltype, CPnt rcen, const vector<CPnt> &cens, 
						const vector<double> &radii, 
						const vector<double> &chg, const vector<CPnt> &cpos, double idiel, 
						const vector<REAL*> &iMats, REAL intraRcutoff,
						const vector<vector<CPnt> > &SPxes,  const  vector<int> &nSPx, 
						const vector<vector<int> >&neighs,
						const vector< vector<int> > &intraPolLists_near,
						const vector<CMulExpan*> &Fself, const vector<CMulExpan*> &Hself,
						const vector<CLocalExpan*> & LFs_intraSelf, 
						const vector<CLocalExpan*> &LHs_intraSelf,
						const vector<CLocalExpan*> & LFs_intraSelf_far, 
						const vector<CLocalExpan*> &LHs_intraSelf_far, 
						const vector<vector<REAL> > &qSolvedFself, 
						const vector<vector<REAL> > &qSolvedHself, 
						const vector<double> &totalFself, const vector<double> &totalHself,
						const vector<CMolCell> &molcell = vector<CMolCell>( 0 ));
	
	// for queries
	CMolecule( CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii, 
						const vector<CMulExpan*> Hself, 
						const vector<CMolCell> &molcell = vector<CMolCell>( 0 ));
	
	~CMolecule(  );
	
	
	static void initConstants( REAL kappa, REAL sdiel );
	static void initMutualConstants( const vector<CMolecule*> & mols, 
																	REAL interRCutoff, REAL interactRCutoff, bool bGrad );
	static void resetMolSystem(  );
	static void deleteConstants(  );
	
	static void aggregateMolMultipoles( vector<CMolecule*> & mols );
	static void aggregateMolMultipoles_Conditional( vector<CMolecule*> & mols );
	
	static void computeMolTypeValues( CPnt rcen, const vector<CPnt> &cens, 
																	 const vector<double> &radii,
																	 const vector<double> &chg, const vector<CPnt> &cpos, double idiel,
																	 REAL intraRcutoff,
																	 const vector<vector<CPnt> > &SPxes, const vector<int> &nSPx,
																	 const  vector<vector<int> >&neighs,
																	 const vector<CMulExpan*> &Fself, const vector<CMulExpan*> &Hself,
																	 const vector<vector<REAL> > &qSolvedFself,
																	 const vector<vector<REAL> > &qSolvedHself,
																	 const vector<double> &totalFself, const vector<double> &totalHself,
																	 const vector< vector<int> > &intraPolLists_near,
																	 const vector< vector<int> > &intraPolLists_far,
																	 vector<CLocalExpan*> &LFs_intraSelf, 
																	 vector<CLocalExpan*> &LHs_intraSelf,
																	 vector<CLocalExpan*> &LFs_intraSelf_far, 
																	 vector<CLocalExpan*> &LHs_intraSelf_far );
	
	static bool computeForces( vector<CMolecule*> & mols,vector<CPnt> & force, 
														vector<CPnt> & torque );
	static REAL computeTotalIntEnergy( const vector<CMolecule*> & mols );
	static REAL computeMolPairPot( vector<CMolecule*> & mols, int i, int j );
	static REAL computeMolIntEnergy( vector<CMolecule*> & mols, int i );
	static REAL computePotInSpace(  CMolecule  mol, CPnt P );
	static REAL computePotInSpace( const vector<CMolecule*> & mols, CPnt P );
	
	static void generateAll( const vector<CPnt> &scens, const vector<double> &srad,
													vector<vector<CPnt> > &SPxes, vector<int> &nSPx,
													vector<vector<int> >&neighs, const vector<CMulExpan*> Fself,
													const vector<CMulExpan*> Hself,
													vector<vector<double> > &qSolvedF,vector<vector<double> > &qSolvedH,
													vector<double> &totalF, vector<double> &totalH,    
													CPnt rcen, vector<CMolCell> &molcells,
													REAL intraRcutoff,
													vector< vector<int> > &intraPolLists_near,
													vector< vector<int> > &intraPolLists_far,
													const vector<double> &chg, const vector<CPnt> &cpos, double idiel,
													vector<CLocalExpan*> &LFs_intraSelf, 
													vector<CLocalExpan*> &LHs_intraSelf,
													vector<CLocalExpan*> &LFs_intraSelf_far,
													vector<CLocalExpan*> &LHs_intraSelf_far );
	
	static void generateInterMolPolList( vector<CMolecule*> & mols );
	static void generateIntraMolPolList( vector<CMolecule*> & mols );
	
	static void generateInterXFormsForPolarize( vector<CMolecule*> & mols );
	static bool generateInterXFormsForPolarize_LowMemory( vector<CMolecule*> & mols );

	//!  The CMolecule generateMolTypeIntraPolLists function
	/*!	This function will generate lists for intra-molecular polarization
	 Add a comment to this line if the distance between sphere ki and kj is 
	 less than intraRcutoff, they will be put into intraPolLists_near, 
	 where the local expansion will be performed numerically ( see Eq. 
	 27a, 28b in 2010 JCTC paper ). Otherwise they will be  put into 
	 intraPolLists_far, where the local expansion will be performed analytically
			\param cens a vector of XYZ coordinates of CG spheres
			\param radii a vector of floating point radii of CG spheres
			\param intraRcutoff a floating point cutoff of the spheres to be considered
							near
			\param intraPolLists_near a vector of spheres that are considered close
			\param intraPolLists_far a vector of spheres that are considered far */	
	static void generateMolTypeIntraPolLists( const vector<CPnt> &cens, 
																					 const vector<double> &radii,
																					 REAL intraRcutoff,
																					 vector< vector<int> > &intraPolLists_near,
																					 vector< vector<int> > &intraPolLists_far );
	
	static void polarize_mutual( vector<CMolecule*> & mols, 
															bool bPot, int farFieldFreq );
	static void prepareDTA( const vector<CMolecule*> & mols, 
												 int j, vector<CGradExpan> &tG );
	static void prepareDTAintra_LowMemory( const vector<CMolecule*> & mols, 
																				int j, vector<CGradExpan> &tG );
	static void prepareDTA_iself_LowMemory( const vector<CMolecule*> & mols, 
																				 int i, int ki );
	static void prepareDTA_LowMemory( const vector<CMolecule*> & mols, 
																	 int j, vector<CGradExpan> &tG );
	
	double recomputeFromSelfVal_LowMemory( const vector<CMolecule*> & mols, 
																				int i, int ki, bool bUpdateFarField );
	double recompute_LowMemory( const vector<CMolecule*> & mols, int i, 
														 int ki, bool bUpdateFarField );
	double recompute_LowMemory( int ki, bool bUpdateFarField );
	
	double recomputeGrad_LowMemory( const vector<CMolecule*> & mols, 
																 const CGradExpan &tG_DTA, 
																 int i, int ki, int j, bool bUpdateFarField,vector<CGradExpan> &tGFs,
																 vector<CGradExpan> &tGHs, vector<CGradExpan> &tGHrots,
																 map<CFullSphereID, int, classCompareCFullSphereID> &gHmap );
	
	void reexpand_LowMemory( const vector<CMolecule*> & mols, 
													int i, int ki, bool bUpdateFarField );
	void reexpand_LowMemory( int ki, bool bUpdateFarField );
	void reexpandLSFromList_LowMemory( const vector<CMolecule*>& mols, int i, int ki );
	void reexpandIntra_Near_LowMemory( int ki, const CLocalExpan &LH0 );
	
	void reexpandFromSelfVal_LowMemory( const vector<CMolecule*> & mols, 
																		 int i, int ki, bool bUpdateFarField );
	void reexpandLSFromSelfVal_LowMemory( const vector<CMolecule*> & mols, 
																			 int i, int ki );
	void reexpandLSFromSelfVal_debug( const vector<CMolecule*> & mols, 
																	 int i, int ki, int j, int kj );

	//!  The CMolecule findNeighbors function
	/*!	Function to determine neighboring CG spheres for sphere ki
			\param ki an integer of the sphere of interest
			\param cens a vector of XYZ coordinates of CG spheres
			\param radii a vector of floating point radii of CG spheres
			\param neighs a vector of neighbor lists */
	static void findNeighbors( int ki, const vector<CPnt> &cens,
														const vector<double> &radii, 
														vector<int> &neigh );
	
	static void generateMolCells( const vector<CPnt>&scens,const vector<double> &srad, 
															 CPnt rcen, vector<CMolCell> &molcells );
															 
	static void generateMolExposedChargesFH( const vector<CMulExpan*> Fself, 
																					const vector<CMulExpan*> Hself, 
																					const vector<CPnt> &scens,const vector<double> &srad, 
																					const vector<vector<CPnt> > &SPEx, const vector<int> & nSPx, 
																					vector<vector<double> > &qSolvedF,vector<vector<double> > &qSolvedH, 
																					vector<double> &totalF, vector<double> &totalH );															 

	
	//!  The CMolecule generateMolSPX function
	/*!	Function to generate xform points for all spheres
			\param scens a vector of XYZ coordinates of the centers of all CG 
								spheres within a mol
			\param srad a vector of floating point radii of CG spheres
			\param SPxes a vector of transforms for the CG spheres 
			\param nSPx a vector of indices for the transforms
			\param neighs a vector of neighbor lists */
	static void generateMolSPX( const vector<CPnt> &scens,const vector<double> &srad, 
														 vector<vector<CPnt> > &SPxes, vector<int> &nSPx,
														 vector<vector<int> >&neighs );	
	//!  The CMolecule getSpherePoints function
	/*!	Function to generate surface points for all spheres, classified
	into buried and exposed
			\param cens a vector of XYZ coordinates of the centers of all CG 
								spheres within a mol
			\param radii a vector of floating point radii of CG spheres
			\param rad2 a vector of radii^2
			\param SP a vector of transforms for the CG sphere (unused?)
			\param NP a vector of indices for the transforms 
			\param ki an integer of the CG sphere of interest
			\param SPE a vector of exposed sphere points for the CG sphere ki
			\param SPB a vector of buried sphere points for the CG sphere ki
			\param neighs a vector of neighbors for the kith sphere */
	static void getSpherePoints( const vector<CPnt> &cens, const vector<double> &radii, 
															const vector<double> &rad2,  const vector<CPnt> &SP, 
															const vector<CPnt> &NP, 
															int ki, vector<CPnt> &SPE,vector<CPnt> &SPB, vector<int> neigh );
	//!  The CMolecule getXFormSpherePoints function
	/*!	Function to generate exposed xform points for all spheres
			\param cens a vector of XYZ coordinates of the centers of all CG 
								spheres within a mol
			\param radii a vector of floating point radii of CG spheres
			\param rad2 a vector of radii^2
			\param SP a vector of transforms for the CG sphere (unused?)
			\param NP a vector of indices for the transforms
			\param neighs a vector of neighbors for the kith sphere 
			\param ki an integer of the CG sphere of interest
			\param SPx a vector of transforms for the CG sphere ki
			\param nSPx a vector of indices for the transforms of kith sphere
			\param Nin a maximum limit for the number of surface points */															
	static void getXFormSpherePoints( const vector<CPnt> &cens, 
																	 const vector<double> &radii, 
																	 const vector<double> &rad2,  
																	 const vector<CPnt> &SP, const vector<CPnt> &NP, 
																	 const vector<int> &neigh, int ki,
																	 vector<CPnt> &SPx, int &nSPx, int Nin );
	
	static bool IsInsideSurface( CPnt P, const vector<CPnt> &SP, 
															const vector<CPnt> &NP );

	
	static REAL maxDist( const CPnt & pnt, const vector<CPnt> & pts );

	static bool OutsideAllBoundaries( const vector<CMolecule*> & mols, CPnt P );

	void reexpandGrad_LowMemory( const vector<CMolecule*> & mols, 
															const CGradExpan &tG_DTA, 
															CGradExpan &tG, int i, int ki, int j, bool bUpdateFarField, 
															const vector<CGradExpan> &tGHrots,
															map<CFullSphereID, int, classCompareCFullSphereID> &gHmap );	
	void reexpandIntraGrad_LowMemory( const vector<CMolecule*> & mols, 
																	 CGradExpan &tGF, CGradExpan &tGH, 
																	 int i, int ki, int j,  bool bUpdateFarField, 
																	 const vector<CGradExpan> &tGFs,
																	 const vector<CGradExpan> &tGHs,
																	 map<CFullSphereID, int, classCompareCFullSphereID> &gHmap  );

	void resetCenters( const vector<CMulExpan> &iF, const vector<CMulExpan> &iH  );

	static void saveConfig( const char * fname, const vector<CMolecule*> & mols );

	
	static void writeMolsPQR( const char * fname, const vector<CMolecule*> & mols );
	static void writeMolsXYZR( const char * fname, const vector<CMolecule*> & mols );
	


	
	void computeMolExposedSurfaceChargesH(  );
	void computeMolExposedSurfaceChargesF(  );
	void computeMolTQ(  );

	//!  The CMolecule computeMaxRadius function
	/*!	Function to compute the maximum radius. Found in moldynamics.cpp
			\param cens a vector of XYZ coordinates of the centers of all CG 
								spheres within a mol
			\param radii a vector of floating point radii of CG spheres
			\return a floating point of the largest radius in collection of spheres  */	
	static double computeMaxRadius( const vector<CPnt> &cens, 
																 const vector<REAL> &radii );
	//!  The CMolecule computeMaxRadius function
	/*!	Function to compute the maximum radius of point charges, plus
	a designated sphere tolerance and set it to the object member maxR 
			\param cpos a vector of XYZ coordinates of the point charges   */	
	void computeMaxRadius( const vector<CPnt> &cpos );
	
	void writeMolExpansions( char* runname ) const;
	
	void computeMol_Force_and_Torque( CPnt &force, CPnt &torque ) const ;
	const REAL computePot(  ) const;
	
	const CQuat getOrient(  ) const {return m_orient;}
	void setPos( const CPnt & newPos ) {  m_rcen = CSystem::pbcPos(newPos); }
	void setOrient( const CQuat &Q ) { m_orient = Q; rotate(CQuat());}
	void translate( const CPnt & trans ) 
	{ m_undorcen = m_rcen; m_rcen = CSystem::pbcPos(m_rcen+trans); }
	void untranslate(  ) { m_rcen = m_undorcen;}
	void rotate( const CQuat & dQ );
	void rotateRotCoeff(  );
	
	bool IsInsideSurface( CPnt P ) const;
	static REAL getC2CDist( const CMolecule *mol1, const CMolecule *mol2 );
	static REAL getC2CDist2( const CMolecule *mol1, const CMolecule *mol2 );
	static REAL getSepDist( const CMolecule *mol1, const CMolecule *mol2 );
	static REAL getS2SDist( const CMolecule *mol1, int ki, 
												 const CMolecule *mol2, int kj );
	static REAL getSepDist( const CMolecule *mol1, int ki, 
												 const CMolecule *mol2, int kj );
	
	static void printMolConfig( const CMolecule *mol, char*fname );
	static bool checkGradSpheres( const vector<CMolecule*> & mols, 
															 int j, int m, int km );
	bool willRotCollide( const vector<CMolecule*> &mols, CQuat dQ ) const ;
	bool willRotCollide_cell( const vector<CMolecule*> &mols, CQuat dQ ) const ;
	bool isCollided( const vector<CMolecule*> &mols ) const ;
	bool isCollided_cell( const vector<CMolecule*> &mols ) const ;
	
	//docking
	//  static void initMaxOverlap( const char* fname );
	
	static void initMolContactRigid( vector<CMolecule*> &mols, int mol1type, 
																	vector<CMolContactRigid> &molcontact );
	static void readMolContactRigid( const char* fname, int &mol1type, 
																	vector<CMolContactRigid> &molcontact, double addDist );
	
	static bool checkDocked( CMolecule* mol1,  CMolecule* molj, bool bDebug  );
	static bool checkDocked_debug( CMolecule* mol1,  CMolecule* molj, bool bDebug  );
	static bool updateDockStatus( vector<CMolecule*> &mols, bool bDebug, int n );
	static void updateDockStatistics( const vector<CMolecule*> &mols, 
																	 char* runname, int n, double t,
																	 vector<ofstream*> &douts, int nSpeciesPerBin );
	static void updateSpecies( const vector<CMolecule*> &mols, 
														int j, vector<bool> &bCounted, int & size );

	void polarize_self( bool bPot, int farFieldFreq );
	void polarize_self( bool bPot ) { polarize_self(bPot, 1); }
	
	/////////////////////////////
	// Public inline functions to print out/set variables!!
	void addDockedNeigh( int j ) {m_dockedNeigh.push_back(j);}

	const bool getbBuried(  ) const {return m_bBuried;}
	const bool getBDocked( int ii ) const {return m_bDocked[ii];}	
	const vector<CPnt> & getCellCens(  ) const {return m_cellCens;}
	REAL getChgSum( ) { return m_chg_sum; }	
	const vector<int> & getDockedNeigh(  ) const {return m_dockedNeigh;}
	const int  getDockedNeighNum(  ) const {return m_dockedNeigh.size();}
	const int getID(  ) const {return m_id;}
	const bool getbInterXForm(  ) const {return m_bInterXForm;}
	const map<CSphereIDPair, int, classCompareCSphereIDPair> & getIntraMap(  ) const {return m_intramap;}
	const REAL getIntraRCutoff(  ) const {return m_intraRcutoff;}
	const CSolExpCenter& getKS( int ki ) const { return *(m_k[ki]);}
	CSolExpCenter& getKS( int ki ) { return *(m_k[ki]);}
	const CSolExpCenter* getpKS( int ki ) const { return m_k[ki]; } 
	CSolExpCenter* getpKS( int ki ) { return m_k[ki]; } 
	const REAL getMaxR(  ) const {return m_maxR;}
	const vector<CMolCell>  getMolCells(  ) const {return m_molcells;}
	const CMolContactRigid &getMolContactRigid( int ii ) const 	
	{return m_molcontactlist[ii];}
	const REAL getMolDev(  ) const {return m_moldev;}
	const CMulExpan& getMolM(  ) const {return m_molM;}
	const REAL getMolTQ(  ) const {return m_molTQ;}
	const int getMolType(  ) const {return m_moltype;}
	const int getNDockSides(  ) const {return m_molcontactlist.size();}
	const int getNKS(  ) const {return m_nks;}			// number of spheres in a molecule (S. Liu)
	const int getOrder(  ) const {return m_p;}	// number of poles (S. Liu)
	const CPnt getRCen(  ) const {return m_rcen;}	
	const int getReciprocalSide( int ii ) const { return getNDockSides()-1 - ii; }
	const CRotCoeff & getRot(  ) const {return m_rot;}
	
	const bool isAggregateM(  ) const {return m_bAggregateM;}

	void resetDockStatus(  );

	void setAggregateM( bool bAgg ) {m_bAggregateM = bAgg;}
	void setBDocked( int ii, bool b ) { m_bDocked[ii] = b;}
	static void setInterRcutoff( REAL rcut ) { m_interRcutoff = rcut; }
	static void setInteractRcutoff( REAL rcut ) { m_interactRcutoff = rcut; }
	void setInterXForm( bool bXFS ) {m_bInterXForm = bXFS;}
	void setMolContactList(  vector<CMolContactRigid> &mlist ) {m_molcontactlist=mlist;}
	void setMolDev( REAL dev  ) { m_moldev = dev;}

	static double m_sdiel;
	static double KAPPA;
	static double m_interRcutoff;
	static double m_interactRcutoff;
	static double INTRA_RCUTOFF;
	static bool m_bInfinite;
	static bool m_bGrad;											//!< A boolean for whether or not to compute the gradient
	static int m_unit;
	static int NUM_POINTS_I;									//!< number of surface points for calculation of Imat
	static int NUM_POINTS_X;									//!< number of surface points for calculation of Imat
	
protected:
	static int N_MOL;													//!< Number of molecules in the system
	static int m_nInterXFS;
	static REAL SPHERETOL;										//!< Sphere tolerance for extracting charges
	static REAL MAXDEV;
	
	CPnt m_rcen;															//!< XYZ coordinates for the molecule's center of geom
	CPnt m_undorcen;													//!< A saved XYZ coord for the molecule's center incase of undo
	double m_idiel;														//!< A floating point of the molecule's dielectric constant
	int m_p;																	//!< An integer of the number of poles of the MPE
	int m_id;																	//!< The integer ID of the molecule for many in a system
	int m_moltype;														//!< The integer ID of the molecule type
	int m_nks;																//!< The number of CG spheres in the molecule
	vector<REAL> m_chg;												//!< A vector of floating point charges in molecule
	double m_chg_sum;													//!< A floating point of the sum of the point charges
	vector<CPnt> m_cpos;											//!< A vector of XYZ coords, one for each point charge
	vector<CPnt> m_atomP;											//!< 
	vector<CPnt> m_SP;												//!< 
	vector<CPnt> m_NP;												//!< 
	vector<int> m_clabels;										//!< 
	REAL m_maxR;															//!< The radius of the largest CG sphere in the molecule
	REAL m_intraRcutoff;											//!< 
	CRotCoeff m_rot;													//!< overall rotation of molecule
	CQuat m_orient;														//!< 
	bool m_bKappa;														//!< 
	bool m_bBuried;														//!< 
	bool m_bAggregateM;												//!< 
	bool m_bInterXForm;												//!< 
	bool m_bOwnIntraXForm;										//!< 
	vector<CFullSphereID> m_interMolPolList;	//!< 
	vector<int> m_intraMolPolList;						//!< 
	vector<int> m_intraMolInteractOnlyList;		//!< 
	
	map<CSphereIDPair, int, classCompareCSphereIDPair> m_intramap;
	
	vector<CSolExpCenter*> m_k;		//!< this contains information about F&H exps, LF, LHN re-expansions( S. Liu )
	CMulExpan m_molM;													//!<
	REAL m_molTQ;															//!<
	REAL m_moldev;														//!<
	
	// cells for clash check
	//const vector<CMolCell> & m_molcells;
	vector<CPnt> m_cellCens;
	
	// for docking purpose
	vector<CMolContactRigid> m_molcontactlist;
	vector<bool> m_bDocked;
	vector<int> m_dockedNeigh;
	
private:

	//! CMolecule assignCharges function
	/* Function to assign charges to respective CG centers
			\param cens a vector of XYZ coordinates of CG spheres
			\param radii a vector of floating point radii of CG spheres
			\param cpos a vector of XYZ coordinates of point charges in the molecule
			\param clabel a vector of integers identifying the numeric ID that each
							point charge is assigned to */	
	static  void assignCharges( const vector<CPnt> &cens, 
														 const vector<double> &radii, 
														 const vector<CPnt> &cpos, vector<int> &clabel );
	void computeMolMultipole(  );														 
	//! CMolecule extractCharges function
	/* Function to extract charge and positions 
			\param ki an integer of the CG sphere of interest
			\param cenKi a vector of XYZ coordinates of the sphere of interest
			\param cpos a vector of XYZ coordinates of point charges in the molecule
			\param clabel a vector of integers identifying the numeric ID that each
							point charge is assigned to
			\param radii a vector of floating point radii of CG spheres
			\param posAssigned a vector of XYZ coordinates of assigned charges WRT sphere ki
			\param chgAssigned a vector of point charges assigned to sphere ki
			\param allPosKi a vector of XYZ of all point charges WRT sphere ki */	
	static void extractCharges( int ki, CPnt cenKi, const vector<CPnt> &cpos, 
														 const vector<double> &chg, const vector<int> & clabel, 
														 vector<CPnt> &posAssigned, vector<double> & chgAssigned );
	//! CMolecule extractCharges function
	/* Function to extract charge and positions and move the point charge to a frame 
	WRT to the charge center that it is assigned
			\param ki an integer of the CG sphere of interest
			\param cenKi a vector of XYZ coordinates of the sphere of interest
			\param cpos a vector of XYZ coordinates of point charges in the molecule
			\param clabel a vector of integers identifying the numeric ID that each
							point charge is assigned to
			\param radii a vector of floating point radii of CG spheres
			\param posAssigned a vector of XYZ coordinates of assigned charges WRT sphere ki
			\param chgAssigned a vector of point charges assigned to sphere ki
			\param allPosKi a vector of XYZ of all point charges WRT sphere ki */		
	static void extractCharges( int ki, CPnt cenKi, const vector<CPnt> &cpos,
														 const vector<double> &chg,const vector<int> & clabel, 
														 vector<CPnt> &posAssigned, vector<double> & chgAssigned, 
														 vector<CPnt> &allPosKi );
														 
														 
	static REAL interactCenters( CMolecule* moli, CMolecule* molj );
	static REAL interactCenters_LowMemory( CMolecule* moli, CMolecule* molj );
	static REAL interactMols( CMolecule* moli, CMolecule* molj );
	static REAL interactMolWithCenters( CMolecule* moli, CMolecule* molj );
	static REAL interactMolWithCenters_LowMemory( CMolecule* moli, 
																							 CMolecule* molj );																							 

	//! CMolecule spherePts function
	/* Function to generate N points uniformly on the surface of a unit sphere
	\param N is the int number of points to generate
	\param th a vector of floating points of theta angles in rads
	\param ph a vector of floating points of phi angles in rads */
	static void spherePts( int N, vector<REAL> & th, vector<REAL> & ph );

	// Printing or inline functions ///////////////////////////////////////////////
	void addInterMolPolList( CFullSphereID id ) {m_interMolPolList.push_back(id);}
	void addIntraMolPolList( int id ) {m_intraMolPolList.push_back(id);}
	void addIntraMolInteractOnlyList( int id ) {m_intraMolInteractOnlyList.push_back(id);}
	
	void clearInterMolPolList(  ) {m_interMolPolList.clear();}
	void clearIntraMolPolLists(  ) {m_intraMolPolList.clear(); m_intraMolInteractOnlyList.clear();}

	const vector<int> & getIntraMolInteractOnlyList() const { return m_intraMolInteractOnlyList;}
	const vector<int> & getIntraMolPolList(  ) const { return m_intraMolPolList;}
	const vector<CFullSphereID> & getInterMolPolList( ) const { return m_interMolPolList;}
	int getInterMolPolListSize(  ) const { return m_interMolPolList.size();}

	//! CMolecule useXFormN function
	/* Function to determine whether or not XFormN for non-overlapping spheres
	use for intra xform only!!!!! not inter because SPX is not rotated 
	\param rho is the floating point distance between the spheres
	\return a bool of whether or not rho is less than the minimum given
						separation distance.  */
	static bool useXFormN( REAL rho ){	return ( rho < MINSEPDIST  );}
		
	// Private variables
	static double m_total;
	
};	// end CMolecule


//inline bool useXFormN( REAL sepdist ){	return ( sepdist < MINSEPDIST  );}	

#endif

