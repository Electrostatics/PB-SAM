#ifndef _CONTACT_H_
#define _CONTACT_H_
#include <vector>
#include "util.h"

using namespace std;

//!  The CContact class
/*! The class that contains information about atomic molecule
contacts */
class CContact
{
public:
	//!  The CContact class constructor
	/*! The constructor that creates a CContact object */
	CContact(int id1, int id2, REAL dist=SEPDIST) : 
	m_id1(id1),m_id2(id2),m_dist(dist) {}

	// functions to print various part of class object
  int getID1() const {return m_id1;}
  int getID2() const {return m_id2;}
  double getDist() const {return m_dist;}
	
	static REAL SEPDIST;					//!< A separation distance
private:
  int m_id1;										//!< The ID of atom on molecule 1 in contact
	int m_id2;										//!< The ID of atom on molecule 1 in contact
  REAL m_dist;									//!< The physical distance between the two atoms when in contact
};

/*#########################################################*/
/*#########################################################*/
// defines the list of contact molecule type 
// makes with another molecule type
/*#########################################################*/
/*#########################################################*/

class CMolContact
{
	public:
		CMolContact( char* fname ); 
		CMolContact( int mol2type, int ncontact, const vector<CContact> &clist );
		CMolContact operator=( const CMolContact &m )
		{
			m_mol2type = m.getMol2Type(  );
			m_ncontact = m.getNContact(  );
			m_bDocked = m.getBDocked(  );
			m_contactlist = m.getContactList(  ); 
		} 

		//CMolContact( char* fname );
		int getMol2Type(  ) const {return m_mol2type;}
		int getNContact(  ) const {return m_ncontact;}
		bool getBDocked(  ) const {return m_bDocked;}
		const vector<CContact> & getContactList(  ) const {return m_contactlist;}

		static void readMolContact( const char* fname, int &mol1type, 
						vector<CMolContact> &molcontact, double sepdist=CContact::SEPDIST );
		static void readMolContactWithDist( const char* fname, int &mol1type, 
						vector<CMolContact> &molcontact  );
		static bool IsDocked( CMolecule* mol1, CMolecule* mol2,
						vector<vector<CMolContact> > MOLCONTACTLIST );

		// static int MINNCONTACT;

		private:
			int m_mol2type;
			int m_ncontact;
			bool m_bDocked;
			vector<CContact> m_contactlist;
};


#endif
