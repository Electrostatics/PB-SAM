#ifndef _HASH_H_
#define _HASH_H_
#include <map>
#include <iostream>

class CSphereID
{
 public:
  CSphereID(int is_) : is(is_) {};
  const int kid() {return is;} 
  int is;
};

////////////////////////////////////////////
class CSphereIDPair
{
 public:
  CSphereIDPair(int is_, int js_) : is(is_), js(js_){};
  int is, js;
};

struct classCompareCSphereIDPair
{ 
  bool operator() (const CSphereIDPair & lhs, const CSphereIDPair &rhs)
  { 
    if(lhs.js <  rhs.js) return true;
    else if (lhs.js >  rhs.js) return false; 
    else //same js -> next compare im
      {
	return (lhs.is <  rhs.is);
      } 
  }//end-operator
};

//////////////////////////////////////////////
class CFullSphereID
{
 public:
  CFullSphereID(int im_, int is_) : im(im_), is(is_) {};
  int mid() const {return im;}
  int kid() const {return is;}
  int im, is;

};

//////////////////////////////////////////////
class CFullSphereIDPair
{
 public:
  CFullSphereIDPair(int im_, int is_, int jm_, int js_) : im(im_), is(is_), jm(jm_), js(js_){};

  int im, is, jm, js;
};

struct classCompareCFullSphereIDPair
{ 
  bool operator() (const CFullSphereIDPair & lhs, const CFullSphereIDPair &rhs)
  { 
    if(lhs.jm <  rhs.jm) return true;
    else if (lhs.jm >  rhs.jm) return false; 
    else //same jm -> next compare js
      {
	if(lhs.js <  rhs.js) return true;
	else if (lhs.js >  rhs.js) return false; 
	else //same js -> next compare im
	  {
	    if(lhs.im <  rhs.im) return true;
	    else if (lhs.im >  rhs.im) return false; 
	    else //same im -> next compare is
	      return (lhs.is <  rhs.is);
	  } 
      }
    
  }//end-operator
};

struct classCompareCFullSphereID
{ 
  bool operator() (const CFullSphereID & lhs, const CFullSphereID &rhs)
  { 
    if(lhs.im <  rhs.im) return true;
    else if (lhs.im >  rhs.im) return false; 
    else //same im -> next compare is
      return (lhs.is <  rhs.is);

  }//end-operator
};
#endif
