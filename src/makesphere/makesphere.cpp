#include <cstdlib>
#include <fstream>
#include <stdio.h>

#include "spheres.h"
#include "readutil.h"

using namespace std;

//!  The mainMakesphere function
/*! The function that uses a MC approach to coarse-graining the input. The function
reads in the PQR and MSMS files.  If this is a contact situation, the function then
reads in the contact atoms and identifies the surface points they correspond to. It then
identifies which spheres are solvent exposed and ensures that they are connected to
a neighboring CG sphere.  Then it prints out the CG spheres and terminates  
		\param pqrfname a path/filename for the input PQR file 
					to be coarse-grained
		\param vertfname the vert file path generated from MSMS
		\param cenfname the centerfile to write out to
		\param tolSP the tolerance of how much the sphere can protrude from SES
		\param makecontacts whether or not to include contacts in 
						CGing. 0 for NO, 1 for YES 
		\param tolIP the interface tolerance only for contact situation 
		\param contactpqrfname the file of contact atom information*/
int mainMakesphere(int argc, char ** argv)
{
  seedRand(-1);
  double start = read_timer();
  const char* pqrfname   = argv[1]; // arg 1 = atompqrfile (in)
  const char* vertfname = argv[2]; // arg 2 = vertfile (in)
  const char* cenfname = argv[3]; // arg 3 = centerfile (out)
  const char* contactpqrfname = argv[7]; // arg 4 = contact.atom.pqrfile (in)
  const double tolSP = atof(argv[4]); // arg 5 = global tolerance (how much sphere can protrude from SES) 
  const double tolIP = atof(argv[6]); // arg 6 = interface tolerance (same as 5, but for interface)

  ofstream fout;
  int np;
   // load PQR files
  vector<CPnt> scen, apos, capos;
  vector<double> srad, CHG, arad, dum, carad;
  readpqrAtoms(pqrfname, apos, CHG ,arad);
  readpqrAtoms(contactpqrfname, capos, dum, carad );
  
  cout <<"Atoms: "<<apos.size()<<endl;  // all atom positions (lab frame)
  // find interface contact positions
  vector <CPnt> conPos;
  const double conCutoff = 5.0;
  getInterfaceAtoms(apos, capos, arad, carad, conPos, conCutoff);
  cout <<"Contact Atoms: "<<conPos.size()<<endl;
  // read molecular surface points from MSMS vertex files
  vector<CPnt> SP, NP;
  readSurface(vertfname, SP, NP);
  // identify interficial vertex points
  vector<int> IPind;
  const double atomVertCutoff = 2.0;
  classifyVertPoints(conPos,SP, IPind, atomVertCutoff );
  // find least number of spheres needed to encompass vertex points
  vector<CPnt> cens;
  vector<double> radii;
 
  findCenters_interface(apos, arad, SP, NP, IPind, cens, radii, tolSP, tolIP);
 // extract solvent centers 
  int nk = cens.size();
  vector<double> rad2(nk);
  for(int ki = 0; ki < nk; ki++) rad2[ki] = radii[ki]*radii[ki] *1.0;
  vector<int> indSCen;
  for(int ki = 0; ki < nk; ki++)
    {
      vector<int> neigh;
      findNeighbors(ki, cens, radii, neigh);
      if( IsSphereExposed(ki, neigh, cens, radii, rad2, SP, NP)  ) indSCen.push_back(ki);
    }
  
  // modify centers to ensure they connect
  cout <<"modifying radius to connect"<<endl;
  bool bAllConnected = false;
  const int maxct = 20;
  const REAL incrementR = 0.5;
  int ct = 0;
  while( !bAllConnected  && ct < maxct )
    {
      bAllConnected = true;
      for(int ki=0; ki<indSCen.size();ki++ )
	{
	  cout <<"ct "<<ct<<" ki "<<ki<<" -> ";
	  vector<int> neigh;
	  findNeighbors(ki, cens, radii, neigh);	
	  if(neigh.size() == 0) 
	    {
	      bAllConnected = false;
	      radii[ki] += incrementR;
	      cout <<" incremented R to "<<radii[ki]<<endl;
	    }
	  else cout <<"connected"<<endl;
	}
      ct++;
    }

  if(ct==maxct) cout <<"warning! reached maxct within completely connected all spheres"<<endl;
  // Write center j and radius out to file
  fout.open(cenfname);  cout<<"Writing centers to "<<cenfname<<endl;
  for(int ki = 0; ki < indSCen.size(); ki++)
    {
      int k = indSCen[ki]; 
      fout<<cens[k].x()<<" "<<cens[k].y()<<" "<<cens[k].z()<<" "<<radii[k]<<endl;
    }
  fout.close();
  double end = read_timer();
  cout <<"total time taken [s] "<<end-start<<endl;
 return 0;   
}


//!  The main function
/*! The function calls the Makesphere routine  
		\param pqrfname a path/filename for the input PQR file 
					to be coarse-grained
		\param vertfname the vert file path generated from MSMS
		\param cenfname the centerfile to write out to
		\param tolSP the tolerance of how much the sphere can protrude from SES
		\param makecontacts whether or not to include contacts in 
						CGing. 0 for NO, 1 for YES 
		\param tolIP the interface tolerance only for contact situation 
		\param contactpqrfname the file of contact atom information*/
int main(int argc, char ** argv)
{
  return mainMakesphere(argc, argv);
}

