#ifndef _CONTACT_H_
#define _CONTACT_H_

#include <vector>
#include <fstream>
#include "util.h"
#include "molecule.h"

using namespace std;

class CContact
{
 public:
  static REAL SEPDIST;

  CContact(int id1, int id2, REAL dist=SEPDIST) : 
    m_id1(id1),m_id2(id2),m_dist(dist) {};

  static void setSepDist(REAL dist) {SEPDIST = dist;} 
  int getID1() const {return m_id1;}
  int getID2() const {return m_id2;}
  REAL getDist() const {return m_dist;}
  //  bool IsContact(vector<CMolecule*> &mols) const;


 private:
  int m_id1, m_id2;
  REAL m_dist;

};


// defines the list of contact molecule type makes with another molecule type
class CMolContact
{

 public:
  CMolContact(char* fname); 
  CMolContact(int mol2type, int ncontact, const vector<CContact> &clist);
  CMolContact operator=(const CMolContact &m)
    {
      m_mol2type = m.getMol2Type();
      m_ncontact = m.getNContact();
      m_bDocked = m.getBDocked();
      m_contactlist = m.getContactList(); 
    } 
    //CMolContact(char* fname);
  int getMol2Type() const {return m_mol2type;}
  int getNContact() const {return m_ncontact;}
  bool getBDocked() const {return m_bDocked;}
  const vector<CContact> & getContactList() const {return m_contactlist;}

  static void readMolContact(const char* fname, int &mol1type, vector<CMolContact> &molcontact, double sepdist=CContact::SEPDIST);
  static void readMolContactWithDist(const char* fname, int &mol1type, vector<CMolContact> &molcontact );
  static bool IsDocked(CMolecule* mol1, CMolecule* mol2,
		       vector<vector<CMolContact> > MOLCONTACTLIST);

  // static int MINNCONTACT;

 private:
  int m_mol2type;
  int m_ncontact;
  bool m_bDocked;
  vector<CContact> m_contactlist;

};


#endif
