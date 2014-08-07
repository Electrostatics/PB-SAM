#ifndef _MOLCELL_H_
#define _MOLCELL_H_

#include <vector>
#include "util.h"

using namespace std;

class CMolCell
{
 public:

  CMolCell() : m_cen(CPnt()),m_rad(0.0) {};
  CMolCell(CPnt cen, REAL rad)  :   m_cen(cen),m_rad(rad) {};
  
  static int getKey(int ix, int iy, int iz) {return ix*4+iy*2+iz;}

  const CPnt & getCen() const {return m_cen;}
  REAL getRad() const {return m_rad;}
  const vector<int> & getSphereList() const {return m_sphereList;}
  int getSphereListSize() const {return m_sphereList.size();}
  void addSphere(int ki) {m_sphereList.push_back(ki);}
  void setRad(double rad) {m_rad = rad; }
 private:
  CPnt m_cen;
  REAL m_rad;
  vector<int> m_sphereList;

};


#endif
