#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include <cfloat>
#include <cassert>

#include "readutil.h"
#include "contact.h"
double CContact::SEPDIST;

using namespace std;

//!  The mainFindContacts_fromAtomContactList function
/*! The function that takes a list of atoms in contact 
and maps to contacting CG spheres
		\param centerpqrfile1 - chain 0 (in) 
		\param centerpqrfile2 - chain 1 (in)
		\param atomcontactlist - expressed as index of apos1 and apos2 (in)
		\param atompqr1fname - chain 0 atom pqr file (in)
		\param atompqr2fname - chain 1 atom pqr file (in)*/
int mainFindContacts_fromAtomContactList(int argc, char ** argv)
{
  cout <<"Chosing centers using atom contact list ..."<<endl;
  seedRand(-1);
  double start = read_timer();
	
	if (argc != 6)
	{
		cout << "Correct usage: ./findContact [cenpqr1fname] "
						<< " [cenpqr2fname] [atomcontactfname] " 
						<< "  [atompqr1fname] [atompqr2fname] " << endl; exit(0);
	}
	
  // load PQR, center files
  const char* cenpqr1fname   = argv[1];
  const char* cenpqr2fname   = argv[2];
  const char* atomcontactfname = argv[3];
  const char* atompqr1fname   = argv[4];
  const char* atompqr2fname   = argv[5];
	
  vector<CPnt> scens1, scens2, apos1, apos2;
  vector<double> srad1, srad2, duma, dumb;
  readpqrAtoms(atompqr1fname, apos1, duma, dumb);
  readpqrAtoms(atompqr2fname, apos2, duma, dumb);
	
  readpqrCenters(cenpqr1fname,scens1,srad1, false);
  readpqrCenters(cenpqr2fname,scens2,srad2, false);
  cout <<"Centers1: "<<scens1.size()<<endl;
  cout <<"Centers2: "<<scens2.size()<<endl;
  
  // read in contact pairs
  vector<CContact> list0;
  ifstream fin(atomcontactfname);
  while (true) 
	{
		int c1, c2;
		double dist;
		fin >> c1 >> c2 ;
		if (fin.eof()) break;
		
		list0.push_back( CContact(c1,c2) );
	} 
	
  cout <<"atom contacts read "<<list0.size()<<endl;
  
  // for each pair of atoms, find which sphere pair they belong to on centerpqrfiles
  // choose the sphere with the smallest radius because this makes 
  // the docking definition more precise	
  vector<CContact> spherelist;
	
  for(int h = 0; h < list0.size(); h++)
	{
		assert(list0[h].getID1() < apos1.size() );
		assert(list0[h].getID2() < apos2.size() );
		
		CPnt a1 = apos1[ list0[h].getID1() ];
		CPnt a2 = apos2[ list0[h].getID2() ];

		int kid1=-1, kid2=-1;
		bool bFound1=false,bFound2=false;
		double minrad = DBL_MAX;
		for(int ki = 0; ki < scens1.size(); ki++)
		{
			double dist2surface = srad1[ki] - (a1-scens1[ki]).norm(); 
			if( dist2surface >= 0)  
	    { 
	      bFound1 = true; 
	      if(srad1[ki] < minrad) { kid1 = ki; minrad = srad1[ki]; }
	    }
		}
		
		minrad = DBL_MAX;
		for(int ki = 0; ki < scens2.size(); ki++)
		{
			double dist2surface = srad2[ki] - (a2-scens2[ki]).norm();
			if( dist2surface >= 0)
			{
				bFound2 = true;
				if(srad2[ki] < minrad) { kid2 = ki; minrad = srad2[ki]; }
			}
		}
		
		spherelist.push_back( CContact(kid1, kid2) );
		
		double dist1 = srad1[kid1] - (a1-scens1[kid1]).norm();
		double dist2 = srad2[kid2] - (a2-scens2[kid2]).norm();
		cout <<kid1<<" "<< dist1 <<" "<<kid2<<" "<< dist2 <<"  " <<dist1+ dist2<<endl;
	}
  
  cout <<"Final Sphere List : "<<spherelist.size()<<endl;
  for(int k=0; k<spherelist.size(); k++)
    cout <<spherelist[k].getID1()<<" "<< spherelist[k].getID2()<<endl;
	
  double end = read_timer();
  cout <<"total time taken [s] "<<end-start<<endl;
  
  return 0;   	
}

//!  The main function
/*! The function that calls mainFindContacts_fromAtomContactList
		\param centerpqrfile1 - chain 0 (in) 
		\param centerpqrfile2 - chain 1 (in)
		\param atomcontactlist - expressed as index of apos1 and apos2 (in)
		\param atompqr1fname - chain 0 atom pqr file (in)
		\param atompqr2fname - chain 1 atom pqr file (in)*/
int main(int argc, char ** argv)
{
	return mainFindContacts_fromAtomContactList(argc, argv);
}
