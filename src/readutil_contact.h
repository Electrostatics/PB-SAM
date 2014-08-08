#include <iostream>
#include <fstream>
#include <vector>
#include "util_contact.h"

// various functions to read files

//!  readpqrAtoms function
/*! A function to read in atoms from a PQR file
 \param fname a character string to read in centers from
 \param pnt a vector of xyz coordinate objects to store coordinates in
 \param ch a vector of charges
 \param R a vector of radii of atoms */
void readpqrAtoms_old(const char * fname, vector<CPnt> & pnt, vector<REAL> & ch,
											vector<REAL> & R)
{
  ifstream fin(fname);
  if (!fin.is_open())
	{
		cout << "Could not open file " << fname << endl;
		exit(0);
	}
	
  pnt.clear();
  ch.clear();
  R.clear();
	
  char buf[100], temp[10];
  int nc = 0;
  REAL sum = 0.0;
	
  fin.getline(buf,99);
	
  while (!fin.eof())
	{
		double x,y,z,c,r;
		if (strncmp(&(buf[0]),"ATOM",4) == 0)
		{
			sscanf(&(buf[31]), "%lf %lf %lf", &x, &y, &z); // position
			sscanf(&(buf[55]), "%lf", &c); // charge
			sscanf(&(buf[63]), "%lf", &r);
			
			pnt.push_back(CPnt(x,y,z));
			ch.push_back(c);
			R.push_back(r);
			sum += c;
			if(fabs(c) > 1e-15) nc++;
		}
		
		fin.getline(buf,99);
	}
  cout << fname << ": Total atoms: " << pnt.size()<<" charged atoms: "<<nc<< ", net charge: "
	<< sum <<endl;
}

//!  readpqrAtoms function
/*! A function to read in atoms from a PQR file
 \param fname a character string to read in centers from
 \param pnt a vector of xyz coordinate objects to store coordinates in
 \param ch a vector of charges
 \param R a vector of radii of atoms */
void readpqrAtoms_new(const char * fname, vector<CPnt> & pnt, vector<REAL> & ch,
											vector<REAL> & R)
{
  ifstream fin(fname);
  if (!fin.is_open())
	{
		cout << "Could not open file " << fname << endl;
		exit(0);
	}
	
  pnt.clear();
  ch.clear();
  R.clear();	
  char buf[100], temp[10];
  int nc = 0;
  REAL sum = 0.0;
	
  fin.getline(buf,99);
	
  while (!fin.eof())
	{
		double x,y,z,c,r;
		if (strncmp(&(buf[0]),"ATOM",4) == 0)
		{
			sscanf(&(buf[31]), "%lf %lf %lf", &x, &y, &z); // position
			sscanf(&(buf[55]), "%lf", &c); // charge
			sscanf(&(buf[62]), "%lf", &r);			
			pnt.push_back(CPnt(x,y,z));
			ch.push_back(c);
			R.push_back(r);
			sum += c;
			if(fabs(c) > 1e-15) nc++;
		}
		
		fin.getline(buf,99);
	}
  cout << fname << ": Total atoms: " << pnt.size()<<" charged atoms: "<<nc<< ", net charge: "
	<< sum <<endl;
}

//!  readpqrCenters function
/*! A function to read in CG spheres from a PQR file
 \param fname a character string to read in centers from
 \param scen a vector of xyz coordinate objects to store coordinates in
 \param srad a vector of floating points of center radii
 \param bNew a boolean of whether or not to use new formatting or not for PQR */
void readpqrCenters(const char * fname, vector<CPnt> & scen, vector<REAL> & srad, bool bNew)
{
  ifstream fin(fname);
  if (!fin.is_open())
	{
		cout << "Could not open file " << fname << endl;
		exit(0);
	}	
  scen.clear();
  srad.clear();	
  int iCoord, iCharge, iRad, iCen;
  if(bNew)
	{
		iCen = 17;
		iCoord = 27;
		iCharge = 52;
		iRad = 59;
	}
  else
	{
		iCen = 18;
		iCoord = 31;
		iCharge = 57;
		iRad = 65;
	}	
  char buf[100], temp[10];	
  fin.getline(buf,99);	
  while (!fin.eof())
	{
		double x,y,z,c,r;
		if (strncmp(&(buf[0]),"ATOM",4) == 0)
		{			
			sscanf(&(buf[iCoord]), "%lf %lf %lf", &x, &y, &z); // position
			sscanf(&(buf[iCharge]), "%lf", &c); // charge
			sscanf(&(buf[iRad]), "%lf", &r);			
			if (strncmp(&(buf[iCen]),"CEN",3) == 0)
	    {
	      srad.push_back(r);
	      scen.push_back(CPnt(x,y,z));
	    }			
		}
		fin.getline(buf,99);
	}
}

