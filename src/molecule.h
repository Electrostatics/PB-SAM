#ifndef _MOLECULE_H_
#define _MOLECULE_H_

#ifdef __OMP
#include "omp.h"
#endif

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
#define MINSEPDIST 5
#define PS_TO_NS 0.001

class CSystem;

class CMolecule
{
 public:

  CMolecule(CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii, const vector<double> &CH, 
	    const vector<CPnt> &POS, double idiel,
	    const vector<vector<CPnt> > &SPxes,  const  vector<int> &nSPx, const vector<vector<int> >&neigh, 
	    const vector< vector<int> > &intraPolLists_near, const vector<CMolCell> &molcell = vector<CMolCell>(0));

  CMolecule(CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii, const vector<double> &CH, 
	    const vector<CPnt> &POS, double idiel,  vector<REAL*> &iMats, REAL intraRcutoff, 
	    const vector<vector<CPnt> > &SPxes,  const  vector<int> &nSPx, const vector<vector<int> >&neighs,
	    const vector< vector<int> > &intraPolLists_near, const vector<CMolCell> &molcell = vector<CMolCell>(0));

  CMolecule(int moltype, CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii, 
	    const vector<double> &chg, const vector<CPnt> &cpos, double idiel, 
	    const vector<REAL*> &iMats, REAL intraRcutoff,
	    const vector<vector<CPnt> > &SPxes,  const  vector<int> &nSPx, const vector<vector<int> >&neighs,
	    const vector< vector<int> > &intraPolLists_near,
	    const vector<CMulExpan*> &Fself, const vector<CMulExpan*> &Hself,
	    const vector<CLocalExpan*> & LFs_intraSelf, const vector<CLocalExpan*> &LHs_intraSelf,
	    const vector<CLocalExpan*> & LFs_intraSelf_far, const vector<CLocalExpan*> &LHs_intraSelf_far, 
	    const vector<vector<REAL> > &qSolvedFself, const vector<vector<REAL> > &qSolvedHself, 
	    const vector<double> &totalFself, const vector<double> &totalHself,
 	    const vector<CMolCell> &molcell = vector<CMolCell>(0));
 
   // for queries
  CMolecule(CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii, 
	    const vector<CMulExpan*> Hself, 
	    const vector<CMolCell> &molcell = vector<CMolCell>(0));

  ~CMolecule();
  // STATIC PUBLIC 
  static void initConstants(REAL kappa, REAL sdiel);
  static void initMutualConstants(const vector<CMolecule*> & mols, REAL interRCutoff, 
				  REAL interactRCutoff, bool bGrad);
  static void resetMolSystem();
  static void deleteConstants();

  static void polarize_mutual(vector<CMolecule*> & mols, bool bPot, int farFieldFreq);
  static void generateInterXFormsForPolarize(vector<CMolecule*> & mols);
  static void generateInterMolPolList(vector<CMolecule*> & mols);
  static void generateIntraMolPolList(vector<CMolecule*> & mols);

  static void aggregateMolMultipoles(vector<CMolecule*> & mols);
  static void aggregateMolMultipoles_Conditional(vector<CMolecule*> & mols);
  static void prepareDTA(const vector<CMolecule*> & mols, int j, vector<CGradExpan> &tG);

  static void generateMolTypeIntraPolLists(const vector<CPnt> &cens, const vector<double> &radii,
					   REAL intraRcutoff,
					   vector< vector<int> > &intraPolLists_near,
					   vector< vector<int> > &intraPolLists_far);
  
  static void computeMolTypeValues(CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii,
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
				   vector<CLocalExpan*> &LFs_intraSelf, vector<CLocalExpan*> &LHs_intraSelf,
				   vector<CLocalExpan*> &LFs_intraSelf_far, vector<CLocalExpan*> &LHs_intraSelf_far);

  static bool generateInterXFormsForPolarize_LowMemory(vector<CMolecule*> & mols);
  static void prepareDTA_LowMemory(const vector<CMolecule*> & mols, int j, vector<CGradExpan> &tG);
  static void prepareDTAintra_LowMemory(const vector<CMolecule*> & mols, int j, vector<CGradExpan> &tG);
  static void prepareDTA_iself_LowMemory(const vector<CMolecule*> & mols, int i, int ki);

  double recomputeFromSelfVal_LowMemory(const vector<CMolecule*> & mols, int i, int ki, bool bUpdateFarField);
  double recompute_LowMemory(const vector<CMolecule*> & mols, int i, int ki, bool bUpdateFarField);
  double recompute_LowMemory(int ki, bool bUpdateFarField);

  void reexpand_LowMemory(const vector<CMolecule*> & mols, int i, int ki, bool bUpdateFarField);
  void reexpand_LowMemory(int ki, bool bUpdateFarField);
  void reexpandLSFromList_LowMemory(const vector<CMolecule*> & mols, int i, int ki);
  void reexpandIntra_Near_LowMemory(int ki, const CLocalExpan &LH0);


  void reexpandFromSelfVal_LowMemory(const vector<CMolecule*> & mols, int i, int ki, bool bUpdateFarField);
  void reexpandLSFromSelfVal_LowMemory(const vector<CMolecule*> & mols, int i, int ki);
  void reexpandLSFromSelfVal_debug(const vector<CMolecule*> & mols, int i, int ki, int j, int kj);

  void reexpandIntraGrad_LowMemory(const vector<CMolecule*> & mols, 
				   CGradExpan &tGF, CGradExpan &tGH, 
				   int i, int ki, int j,  bool bUpdateFarField, 
				   const vector<CGradExpan> &tGFs,
				   const vector<CGradExpan> &tGHs,
				   map<CFullSphereID, int, classCompareCFullSphereID> &gHmap );

  void reexpandGrad_LowMemory(const vector<CMolecule*> & mols, const CGradExpan &tG_DTA, 
		    CGradExpan &tG, int i, int ki, int j, bool bUpdateFarField, 
			      const vector<CGradExpan> &tGHrots,
			      map<CFullSphereID, int, classCompareCFullSphereID> &gHmap);
  double recomputeGrad_LowMemory(const vector<CMolecule*> & mols, const CGradExpan &tG_DTA, 
				 int i, int ki, int j, bool bUpdateFarField,vector<CGradExpan> &tGFs,
				 vector<CGradExpan> &tGHs, vector<CGradExpan> &tGHrots,
				 map<CFullSphereID, int, classCompareCFullSphereID> &gHmap);

  static bool computeForces(vector<CMolecule*> & mols,vector<CPnt> & force, vector<CPnt> & torque);
  static REAL computeTotalIntEnergy(const vector<CMolecule*> & mols);
  static REAL computeMolPairPot(vector<CMolecule*> & mols, int i, int j);
  static REAL computeMolIntEnergy(vector<CMolecule*> & mols, int i);
  static REAL computePotInSpace(const vector<CMolecule*> & mols, CPnt P);
  static REAL maxDist(const CPnt & pnt, const vector<CPnt> & pts);

  static void generateMolCells(const vector<CPnt> &scens,const vector<double> &srad, CPnt rcen, vector<CMolCell> &molcells);
  static void generateMolSPX(const vector<CPnt> &scens,const vector<double> &srad, 
			     vector<vector<CPnt> > &SPxes, vector<int> &nSPx,
			     vector<vector<int> >&neighs);
  static void generateMolExposedChargesFH(const vector<CMulExpan*> Fself, const vector<CMulExpan*> Hself, 
				       const vector<CPnt> &scens,const vector<double> &srad, 
				       const vector<vector<CPnt> > &SPEx, const vector<int> & nSPx, 
				       vector<vector<double> > &qSolvedF,vector<vector<double> > &qSolvedH, 
				       vector<double> &totalF, vector<double> &totalH);
  static void findNeighbors(int ki, const vector<CPnt> &cens,const vector<double> &radii, 
			    vector<int> &neigh);
  static void getSpherePoints(const vector<CPnt> &cens, const vector<double> &radii, 
			      const vector<double> &rad2,  const vector<CPnt> &SP, const vector<CPnt> &NP, 
			      int ki, vector<CPnt> &SPE,vector<CPnt> &SPB, vector<int> neigh);
  static void getXFormSpherePoints(const vector<CPnt> &cens, const vector<double> &radii, 
				     const vector<double> &rad2,  const vector<CPnt> &SP, const vector<CPnt> &NP, 
				     const vector<int> &neigh, int ki,
				     vector<CPnt> &SPx, int &nSPx, int Nin);
  
  static bool IsInsideSurface(CPnt P, const vector<CPnt> &SP, const vector<CPnt> &NP);
  static bool OutsideAllBoundaries(const vector<CMolecule*> & mols, CPnt P);
  static void setInterRcutoff(REAL rcut) { m_interRcutoff = rcut; }
  static void setInteractRcutoff(REAL rcut) { m_interactRcutoff = rcut; }
  static void saveConfig(const char * fname, const vector<CMolecule*> & mols);

  static void writeMolsPQR(const char * fname, const vector<CMolecule*> & mols);

  static double m_sdiel, KAPPA, m_interRcutoff, m_interactRcutoff, INTRA_RCUTOFF;
  static bool m_bInfinite, m_bGrad;
  static int m_unit, NUM_POINTS_I, NUM_POINTS_X;
  //  static vector<REAL> MAXSPHERE_OVERLAP; 
  //void setOrder(const int p);

  void setAggregateM(bool bAgg) {m_bAggregateM = bAgg;}
  void setInterXForm(bool bXFS) {m_bInterXForm = bXFS;}
  //void setIntraRcutoff(REAL rcut) { m_intraRcutoff = rcut; } //debug
  void setMolDev(REAL dev ) { m_moldev = dev;}

  const int getNKS() const {return m_nks;} //number of spheres in a molecule (S. Liu)
  const bool getbBuried() const {return m_bBuried;}
  const bool getbInterXForm() const {return m_bInterXForm;}
  const bool isAggregateM() const {return m_bAggregateM;}
  const CPnt getRCen() const {return m_rcen;}
  const int getOrder() const {return m_p;} //number of poles (S. Liu)

  const int getID() const {return m_id;}
  const int getMolType() const {return m_moltype;}
  const vector<CMolCell>  getMolCells() const {return m_molcells;}
  const vector<CPnt> & getCellCens() const {return m_cellCens;}

  const REAL getMaxR() const {return m_maxR;}
  const REAL getMolDev() const {return m_moldev;}
  const REAL getMolTQ() const {return m_molTQ;}
  const CMulExpan& getMolM() const {return m_molM;}
  const CRotCoeff & getRot() const {return m_rot;}

  void resetCenters(const vector<CMulExpan> &iF, const vector<CMulExpan> &iH );
  const map<CSphereIDPair, int, classCompareCSphereIDPair> & getIntraMap() const {return m_intramap;}

  const REAL getIntraRCutoff() const {return m_intraRcutoff;}

  const CSolExpCenter& getKS(int ki) const { return *(m_k[ki]);}
  CSolExpCenter& getKS(int ki) { return *(m_k[ki]);}
  const CSolExpCenter* getpKS(int ki) const { return m_k[ki]; } 
  CSolExpCenter* getpKS(int ki) { return m_k[ki]; } 

  void polarize_self(bool bPot, int farFieldFreq);
  void polarize_self(bool bPot) { polarize_self(bPot, 1); }

  void computeMolExposedSurfaceChargesH();
  void computeMolExposedSurfaceChargesF();
  void computeMolTQ();

  static double computeMaxRadius(const vector<CPnt> &cens, const vector<REAL> &radii);
  void computeMaxRadius(const vector<CPnt> &cpos);

  void writeMolExpansions(char* runname) const;


  void computeMol_Force_and_Torque(CPnt &force, CPnt &torque) const ;
  const REAL computePot() const;

  const CQuat getOrient() const {return m_orient;}
  void setPos(const CPnt & newPos) { m_rcen = CSystem::pbcPos(newPos); }
  void setOrient(const CQuat &Q) { m_orient = Q; rotate(CQuat());}
  void translate(const CPnt & trans) { m_undorcen = m_rcen; m_rcen = CSystem::pbcPos(m_rcen+trans); }
  void untranslate() { m_rcen = m_undorcen;}
  void rotate(const CQuat & dQ);
  void rotateRotCoeff();

  bool IsInsideSurface(CPnt P) const;
  static REAL getC2CDist(const CMolecule *mol1, const CMolecule *mol2);
  static REAL getC2CDist2(const CMolecule *mol1, const CMolecule *mol2);
  static REAL getSepDist(const CMolecule *mol1, const CMolecule *mol2);
  static REAL getS2SDist(const CMolecule *mol1, int ki, const CMolecule *mol2, int kj);
  static REAL getSepDist(const CMolecule *mol1, int ki, const CMolecule *mol2, int kj);
  
  static void printMolConfig(const CMolecule *mol, char*fname);
  static bool checkGradSpheres(const vector<CMolecule*> & mols, int j, int m, int km);
  bool willRotCollide(const vector<CMolecule*> &mols, CQuat dQ) const ;
  bool willRotCollide_cell(const vector<CMolecule*> &mols, CQuat dQ) const ;
  bool isCollided(const vector<CMolecule*> &mols) const ;
  bool isCollided_cell(const vector<CMolecule*> &mols) const ;
  
  //docking
  //  static void initMaxOverlap(const char* fname);
  
  static void initMolContactRigid(vector<CMolecule*> &mols, int mol1type, vector<CMolContactRigid> &molcontact);
  static void readMolContactRigid(const char* fname, int &mol1type, vector<CMolContactRigid> &molcontact, double addDist);
  
  static bool checkDocked(CMolecule* mol1,  CMolecule* molj, bool bDebug );
  static bool checkDocked_debug(CMolecule* mol1,  CMolecule* molj, bool bDebug );
  static bool updateDockStatus(vector<CMolecule*> &mols, bool bDebug, int n);
  static void updateDockStatistics(const vector<CMolecule*> &mols, char* runname, int n, double t,
				   vector<ofstream*> &douts, int nSpeciesPerBin);
  static void updateSpecies(const vector<CMolecule*> &mols, int j, vector<bool> &bCounted, int & size);

  const vector<int> & getDockedNeigh() const {return m_dockedNeigh;}
  const int  getDockedNeighNum() const {return m_dockedNeigh.size();}
  
  const int getNDockSides() const {return m_molcontactlist.size();}
  const CMolContactRigid &getMolContactRigid(int ii) const {return m_molcontactlist[ii];}
  const int getReciprocalSide(int ii) const { return getNDockSides()-1 - ii; }
  const bool getBDocked(int ii) const {return m_bDocked[ii];}
  void resetDockStatus();
  void setMolContactList( vector<CMolContactRigid> &mlist) {m_molcontactlist=mlist;}
  void setBDocked(int ii, bool b) { m_bDocked[ii] = b;}
  void addDockedNeigh(int j) {m_dockedNeigh.push_back(j);}


 protected:

  static int N_MOL;
  static int m_nInterXFS;
  static REAL SPHERETOL;
  static REAL MAXDEV;

  CPnt m_rcen, m_undorcen;
  double m_idiel;
  int m_p, m_id, m_moltype;
  int m_nks;
  vector<REAL> m_chg;
  vector<CPnt> m_cpos, m_atomP, m_SP, m_NP;
  vector<int> m_clabels;
  REAL m_maxR, m_intraRcutoff;
  CRotCoeff m_rot; // overall rotation of molecule
  CQuat m_orient;
  bool m_bKappa;
  bool m_bBuried;
  bool m_bAggregateM;
  bool m_bInterXForm, m_bOwnIntraXForm;
  vector<CFullSphereID> m_interMolPolList;
  vector<int> m_intraMolPolList, m_intraMolInteractOnlyList;

  map<CSphereIDPair, int, classCompareCSphereIDPair> m_intramap;
    
  vector<CSolExpCenter*> m_k;  //this contains information about F and H expansions, LF, LHN re-expansions(S. Liu)
  CMulExpan m_molM; 
  REAL m_molTQ;
  REAL m_moldev;

  // cells for clash check
  const vector<CMolCell> & m_molcells;
  vector<CPnt> m_cellCens;
  
 // for docking purpose
  vector<CMolContactRigid> m_molcontactlist;
  vector<bool> m_bDocked;
  vector<int> m_dockedNeigh;
  
 private:

  static void spherePts(int N, vector<REAL> & th, vector<REAL> & ph);

  static REAL interactCenters_LowMemory(CMolecule* moli, CMolecule* molj);
  static REAL interactMolWithCenters_LowMemory(CMolecule* moli, CMolecule* molj);
  static REAL interactCenters(CMolecule* moli, CMolecule* molj);
  static REAL interactMolWithCenters(CMolecule* moli, CMolecule* molj);
  static REAL interactMols(CMolecule* moli, CMolecule* molj);
  static bool useXFormN(REAL rho);

  static double m_total;
  
  static  void assignCharges(const vector<CPnt> &cens, const vector<double> &radii, const vector<CPnt> &cpos, 
			     vector<int> &clabel);
  static void extractCharges(int ki, CPnt cenKi, const vector<CPnt> &cpos, 
		      const vector<double> &chg, const vector<int> & clabel, 
		      vector<CPnt> &posAssigned, vector<double> & chgAssigned);
  static void extractCharges(int ki, CPnt cenKi, const vector<CPnt> &cpos,
                      const vector<double> &chg,const vector<int> & clabel, 
		      vector<CPnt> &posAssigned, vector<double> & chgAssigned, vector<CPnt> &allPosKi);
  void clearInterMolPolList() {m_interMolPolList.clear();}
  void addInterMolPolList(CFullSphereID id) {m_interMolPolList.push_back(id);}
  void clearIntraMolPolLists() {m_intraMolPolList.clear(); m_intraMolInteractOnlyList.clear();}
  void addIntraMolPolList(int id) {m_intraMolPolList.push_back(id);}
  void addIntraMolInteractOnlyList(int id) {m_intraMolInteractOnlyList.push_back(id);}
  const vector<CFullSphereID> & getInterMolPolList() const { return m_interMolPolList;}
  const vector<int> & getIntraMolPolList() const { return m_intraMolPolList;}
  const vector<int> & getIntraMolInteractOnlyList() const { return m_intraMolInteractOnlyList;}
  int getInterMolPolListSize() const { return m_interMolPolList.size();}
  void computeMolMultipole();

};
/*
// set the maximum pole and rot accordingly
inline void
CMolecule::setOrder(const int p)
{
  assert (p <= N_POLES);
  assert (m_rot.getOrder() == m_p);

  if (m_p < p)  
    while (m_rot.getOrder() < p) 
      {
	m_rot.incOrder();
	//	cout << "incrementing molecule m_rot: " << m_rot.getOrder()<<endl;
      }
  else if( m_p > p)
    for (int k = m_p; k > p; k--)	m_rot.decOrder(); 

  m_p = p;

  //  assert (m_rot.getOrder() == m_p);
}
*/


// criteria for using CXFormN for non-overlapping spheres
// use for intra xform only!!!!! not inter because SPX is not rotated
inline bool
CMolecule::useXFormN(REAL sepdist)
{
  return (sepdist < MINSEPDIST );
}



#endif
