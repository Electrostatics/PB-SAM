#ifndef _BDNAM_H_
#define _BDNAM_H_
#include <vector>
#include "util.h"
#include "readutil.h"
#include "molecule.h"
#include "contact.h"
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
#define b_DIST 100.0
#define q_DIST 500.0
#define f_DIST 100.0

using namespace std;

// general functions
class CBDnam
{
 public:

  CBDnam(vector<char*> molfnames1, vector<char*> molfnames2, REAL idiel);
  ~CBDnam() {  /*
    for(int i=0; i<m_mols.size(); i++) delete [] m_mols[i];
    for(int i=0; i<m_iMats.size(); i++) delete [] m_iMats[i];
	    */
  }

  enum STATUS {ESCAPED, DOCKED, RUNNING, STUCK};
  static const char STATUS_STRING[4][10];
  static void initConstants(int nMolType);
  static void addMolContact(int mol1type, const vector<CMolContact>&molcontactlist);

  CBDnam::STATUS run(int &scount, char* fname);
  CMolecule* getMol(int i) {return m_mols[i];}
  static CPnt getRandVec(REAL std)
    { return std*CPnt(normRand(), normRand(), normRand()); }


 private:
  // void saveState();
  bool IsDocked(CMolecule* mol1, CMolecule* mol2) const;
  bool escaped(REAL dist) const;
  REAL compute_dt() const; 
  CBDnam::STATUS makeMove(const CPnt & dR2, const CPnt & dO1, const CPnt & dO2, REAL dt);
  

  static int m_nMolType;
  static vector< vector<CMolContact> > MOLCONTACTLIST;

  vector<CMolecule*> m_mols;

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


  // for nam bd studies
  REAL m_Dtr;
  REAL m_Dr1, m_Dr2;
  CQuat m_rot1, m_rot2;


};

inline bool
CBDnam::escaped(REAL dist) const
{
  return (CMolecule::getC2CDist(m_mols[0], m_mols[1]) > dist);
}


#endif
