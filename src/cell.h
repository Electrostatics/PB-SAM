#ifndef _MOLCELL_H_
#define _MOLCELL_H_

#include <vector>
#include "util.h"

using namespace std;

//!  The CMolCell class
/*! The class that contains information about a cell.
 Not sure if it is used yet, but it was meant to speed up
 calculations. */
class CMolCell
{
public:
	//!  The CContact class constructor
	/*! The constructor that creates a CMolCell object
	 and initializes everything at zero. */
  CMolCell() : m_cen(CPnt()),m_rad(0.0) {};
	//!  The CContact class constructor
	/*! The constructor that creates a CMolCell object
	 \param cen a vector of XYZ coords for the cell
	 \param rad a floating point of the radius of the cell. */
  CMolCell(CPnt cen, REAL rad)  :   m_cen(cen),m_rad(rad) {};
  
  static int getKey(int ix, int iy, int iz) {return ix*4+iy*2+iz;}
	
  const CPnt & getCen() const {return m_cen;}
  REAL getRad() const {return m_rad;}
  const vector<int> & getSphereList() const {return m_sphereList;}
  int getSphereListSize() const {return m_sphereList.size();}
  void addSphere(int ki) {m_sphereList.push_back(ki);}
  void setRad(double rad) {m_rad = rad; }
private:
	CPnt m_cen;								//!< XYZ coords of the center of the cell
	REAL m_rad;								//!< Radius of the cell
	vector<int> m_sphereList;	//!< List of spheres in the cell
}; // end CMolCell


#endif
