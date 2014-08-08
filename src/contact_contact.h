#ifndef _CONTACT_H_
#define _CONTACT_H_
#include <vector>
#include "util_contact.h"

using namespace std;

class CContact
{
public:
  static REAL SEPDIST;
  int getID1() const {return m_id1;}
  int getID2() const {return m_id2;}
  double getDist() const {return m_dist;}
  CContact(int id1, int id2, REAL dist=SEPDIST) :
	m_id1(id1),m_id2(id2),m_dist(dist) {}
private:
  int m_id1, m_id2;
  REAL m_dist;
	
};


#endif
