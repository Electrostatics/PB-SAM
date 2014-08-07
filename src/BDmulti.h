#ifndef _BDM_H_
#define _BDM_H_
#include <vector>
#include "util.h"
#include "readutil.h"
#include "molecule.h"
#include "system.h"

#define FACT2 5.604586e2 //(COUL_K*IKbT)
const double AVOGADRO_NUM = 6.0221415e23;
const double ELECT_CHG = 1.60217646e-19;
const double PERMITTIVITY_VAC = 8.854187817e-12;
const double KB = 1.3806503e-23;
const double ANGSTROM = 1e-10;
const double LITRE = 1e-3;
const double IKbT = 1.0/(KB*298.15);
const double COUL_K = 2.30708e-18; // (1/4PI*e0) * e2 / ANGSTROM 

// SYSTEM SPECIFIC PARAMETERS

using namespace std;

// general functions
class CBDmulti
{
 public:

  /*  CBDmulti(int Np, vector<char*> molfnames, char* configfname, REAL idiel,	 
      char* contactfile, double contactSepDist);*/
  CBDmulti(int Np1, int Np2, const vector<char*> &molfnames1, 
	   const vector<char*> &molfnames2, REAL idiel);
  //  CBDmulti(vector<char*> molfnames, REAL idiel, char* contactfile, double contactSepDist); // debug

  ~CBDmulti() {  
    /*    for(int i=0; i<m_mols.size(); i++) delete m_mols[i];
    for(int i=0; i<m_iF.size(); i++) 
    {
	delete m_iF[i];
	delete m_iH[i];
	delete [] m_iMats[i];
      }
    */
  }

  static void initConstants(int nMolType);
  static void addMolContact(int mol1type, const vector<CMolContact>&molcontactlist);



  void resetLattice();
  void restartConfig(char* configname);

  double runFirstDockTime(double maxtime, char* fname);
  CMolecule* getMol(int i) {return m_mols[i];}
  static CPnt getRandVec(REAL std)
    { return std*CPnt(normRand(), normRand(), normRand()); }

  vector<CMolecule*> m_mols; // debug in public
  vector<vector<CPnt> > m_SPxes;
  vector<CMolCell> m_molcells;

 private:
  // void saveState();
  bool IsDocked(CMolecule* mol1, CMolecule* mol2) const;
  
  void initMolTypes(const vector<char*> &molfnames1,
		    const vector<char*> &molfnames2, REAL idiel);
  
  double propagate(double mindist, bool bWrite); // debug
  double propagate(double mindist) { return propagate(mindist, false);}

  REAL computeMinSepDist();
  REAL compute_dt(double minsepdist) const; 
  void makeMove(const vector<CPnt> & dR, const vector<CPnt> & dO, REAL dt);

  //void computeBondingForce(vector<CPnt> &force, vector<CPnt> &torque) const;    
  //void computeBondingForceRigid(vector<CPnt> &force, vector<CPnt> &torque) const;    

  static REAL MINDIST, MAXDIST;
  
  static int m_nMolType;
  static vector< vector<CMolContact> > MOLCONTACTLIST;

  // moltype variables                                                                                
  vector<REAL*> m_iMats1, m_iMats2;
  vector<CMulExpan*> m_iF1,m_iH1, m_iF2,m_iH2;
  vector<vector<double> > m_qSolvedF1, m_qSolvedH1, m_qSolvedF2, m_qSolvedH2;
  vector<double> m_totalF1, m_totalH1, m_totalF2, m_totalH2;
  vector<CLocalExpan*> m_LFs_intraSelf1, m_LHs_intraSelf1, m_LFs_intraSelf_far1, m_LHs_intraSelf_far1;
  vector<CLocalExpan*> m_LFs_intraSelf2, m_LHs_intraSelf2, m_LFs_intraSelf_far2, m_LHs_intraSelf_far2;

  vector<CMolCell> m_molcells1, m_molcells2;

  vector<vector<int> > m_neighs1, m_neighs2;
  vector<vector<CPnt> > m_SPxes1, m_SPxes2;
  vector<int> m_nSPx1, m_nSPx2;
  vector< vector<int> > m_intraPolLists_near1, m_intraPolLists_near2;

  // BD variables
  int m_np1, m_np2;
  vector<REAL> m_Dtr;
  vector<REAL> m_Dr;
  vector<CQuat> m_rots;
  CPnt m_initialrcen1, m_initialrcen2;
};

#endif
