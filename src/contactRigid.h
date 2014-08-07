#ifndef _CONTACTR_H_
#define _CONTACTR_H_

#include <vector>
#include <fstream>
#include "util.h"
using namespace std;

// defines the list of contact molecule type makes with another molecule type
class CMolContactRigid
{

 public:
  CMolContactRigid(int mol2type, double rideal2, CPnt o1A, CPnt o1B, CPnt o2A, CPnt o2B ) : 
    m_mol2type(mol2type), m_rideal2(rideal2),m_o1A(o1A), m_o2A(o2A), m_o1B(o1B), m_o2B(o2B) 
    {
      cout << m_mol2type<<" "<<m_rideal2<<" "
	   <<m_o1A<<" "<<m_o2A<<" "<<m_o1B<<" "<<m_o2B<<endl;
    }
  
  CMolContactRigid operator=(const CMolContactRigid &m)
    {
      m_mol2type = m.getMol2Type();
      m_o1A = m.getVec1A();
      m_o2A = m.getVec2A();  
      m_o1B = m.getVec1B();
      m_o2B = m.getVec2B();
      m_rideal2 = m.getRideal2();
      m_bDocked = m.getBDocked();
    } 
    //CMolContact(char* fname);
  int getMol2Type() const {return m_mol2type;}
  double getRideal2() const {return m_rideal2;}
  bool getBDocked() const {return m_bDocked;}
  const CPnt & getVec1A() const {return m_o1A;}
  const CPnt & getVec1B() const {return m_o1B;}
  const CPnt & getVec2A() const {return m_o2A;}
  const CPnt & getVec2B() const {return m_o2B;}

  static double angleTolerance, maxC2CDist2; 

 private:

  int m_mol2type;
  int m_ncontact;
  bool m_bDocked;
  CPnt m_o1A, m_o2A, m_o1B, m_o2B;
  REAL m_rideal2;

};


#endif
