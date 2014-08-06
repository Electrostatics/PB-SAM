#ifndef _CONTACTR_H_
#define _CONTACTR_H_

#include <vector>
#include <fstream>
#include "util.h"
using namespace std;

//!  The CMolContactRigid class
/*! The class that contains information about the list of
 contact molecule type makes with another molecule type */
class CMolContactRigid
{
public:
	//!  The CMolContactRigid class constructor
	/*! The constructor that creates a CMolContactRigid object
	 \param mol2type the integer ID of atom on molecule 2 in contact
	 \param rideal2 a floating point of the physical distance between the two
	 atoms when in contact 
	 \param o1a a vector of contact details 
	 \param o1b a vector of contact details
	 \param o2a a vector of contact details
	 \param o2b a vector of contact details */
	CMolContactRigid( int mol2type, double rideal2, CPnt o1A,
									 CPnt o1B, CPnt o2A, CPnt o2B  ) :
	m_mol2type( mol2type ), m_rideal2(rideal2), m_o1A(o1A),
	m_o2A( o2A ), m_o1B(o1B), m_o2B(o2B)
	{
		cout << m_mol2type<<" "<<m_rideal2<<" "
		<<m_o1A<<" "<<m_o2A<<" "<<m_o1B<<" "<<m_o2B<<endl;
	}
	
	CMolContactRigid operator=( const CMolContactRigid &m )
	{
		m_mol2type = m.getMol2Type(  );
		m_o1A = m.getVec1A(  );
		m_o2A = m.getVec2A(  );
		m_o1B = m.getVec1B(  );
		m_o2B = m.getVec2B(  );
		m_rideal2 = m.getRideal2(  );
		m_bDocked = m.getBDocked(  );
	}

	// Printing function
	int getMol2Type(  ) const {return m_mol2type;}
	double getRideal2(  ) const {return m_rideal2;}
	bool getBDocked(  ) const {return m_bDocked;}
	const CPnt & getVec1A(  ) const {return m_o1A;}
	const CPnt & getVec1B(  ) const {return m_o1B;}
	const CPnt & getVec2A(  ) const {return m_o2A;}
	const CPnt & getVec2B(  ) const {return m_o2B;}
	
	static double angleTolerance;
	static double maxC2CDist2;
	
private:
	int m_mol2type;	//!< The type of the molecule in contact with
	int m_ncontact;	//!< The number of contact points
	bool m_bDocked;	//!< A boolean of whether it is docked or not
	CPnt m_o1A;			//!< A vector 1A
	CPnt m_o2A;			//!< A vector 2A
	CPnt m_o1B;			//!< A vector 1B
	CPnt m_o2B;			//!< A vector 2B
	REAL m_rideal2;	//!< A real of the distance between mols
}; // end CMolContactRigid

#endif
