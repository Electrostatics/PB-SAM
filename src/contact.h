#ifndef _CONTACT_H_
#define _CONTACT_H_

#include <vector>
#include <fstream>
#include "util.h"
#include "molecule.h"

using namespace std;

//!  The CContact class
/*! The class that contains information about atomic molecule
 contacts */
class CContact
{
public:
  static REAL SEPDIST;					//!< A separation distance
	//!  The CContact class constructor
	/*! The constructor that creates a CContact object
	 \param id1 the integer ID of atom on molecule 1 in contact
	 \param id2 the integer ID of atom on molecule 2 in contact
	 \param dist a floating point of the physical distance between the two
	 atoms when in contact */
  CContact(int id1, int id2, REAL dist=SEPDIST) :
	m_id1(id1),m_id2(id2),m_dist(dist) {};
	
  static void setSepDist(REAL dist) {SEPDIST = dist;}
  int getID1() const {return m_id1;}
  int getID2() const {return m_id2;}
  REAL getDist() const {return m_dist;}
	
private:
  int m_id1;										//!< The ID of atom on molecule 1 in contact
	int m_id2;										//!< The ID of atom on molecule 1 in contact
  REAL m_dist;									//!< The physical distance between the two atoms when in contact
};


//!  The CMolContact class
/*! The class that contains information about molecule
 contacts. It defines the list of contact molecule type makes
 with another molecule type */
class CMolContact
{
public:
	//!  The CMolContact class constructor
	/*! The constructor that creates a CMolContact object
	 \param fname a character string of an input contact filename */
  CMolContact(char* fname);
	//!  The CMolContact class constructor
	/*! The constructor that creates a CMolContact object
	 \param mol2type an integer of the ID of molecule type to which it
	 will be contacted
	 \param ncontact an integer of the number of contact sites in the molecule
	 \param clist a vector of CContact object, one for each pair of atoms in
	 contact */
  CMolContact(int mol2type, int ncontact, const vector<CContact> &clist);
	//!  The CMolContact = operator
	/*! The constructor that creates a CMolContact object
	 from another. */
  CMolContact operator=(const CMolContact &m)
	{
		m_mol2type = m.getMol2Type();
		m_ncontact = m.getNContact();
		m_bDocked = m.getBDocked();
		m_contactlist = m.getContactList();
	}
  int getMol2Type() const {return m_mol2type;}
  int getNContact() const {return m_ncontact;}
  bool getBDocked() const {return m_bDocked;}
  const vector<CContact> & getContactList() const {return m_contactlist;}
	//!  The CMolContact readMolContact function
	/*! The function that reads in the molecular contacts
	 \param fname a character string of a contact file
	 \param mol1type an ID of the first molecule for contact
	 \param molcontact a vector of molecular contacts for the first molecule
	 \param sepdist a floating point of distance in A for separation criteria */
  static void readMolContact(const char* fname, int &mol1type, vector<CMolContact> &molcontact, double sepdist=CContact::SEPDIST);
	//!  The CMolContact readMolContactWithDist function
	/*! The function that reads in the molecular contacts as
	 well as their separation criteria
	 \param fname a character string of a contact file
	 \param mol1type an ID of the first molecule for contact
	 \param molcontact a vector of molecular contacts for the first molecule */
  static void readMolContactWithDist(const char* fname, int &mol1type, vector<CMolContact> &molcontact );
	//!  The CMolContact IsDocked function
	/*! The function that determines whether or not two molecules have docked
	 \param mol1 a CMolecule for the first mol that may be in contact
	 \param mol2 a CMolecule for the 2nd mol that may be in contact
	 \param MOLCONTACTLIST an array of docking criteria for each pair of molecules */
  static bool IsDocked(CMolecule* mol1, CMolecule* mol2,
											 vector<vector<CMolContact> > MOLCONTACTLIST);
	
private:
	int m_mol2type;									//!< The numerical type of the second molecule in the contact pair
	int m_ncontact;									//!< The number of contacts (atom pairs) in the molecule docking
	bool m_bDocked;									//!< A boolean to indicate whether the molecule is docked or not
	vector<CContact> m_contactlist; //!< A list of CContact objects that contain the list of atom pairs in contact
};


#endif
