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
int mainMakesphere( int argc, char ** argv )
{
	// Random number seed
	seedRand( -1 );
	// Timer
	double start = read_timer(  );
	
	if ( (argc != 6) && (argc != 8 ) )
	{
		cout << "Correct usage: ./makesphere [pqrfname] "
				<< " [vertfname] [cenfname] [tolSP] [makecontacts] [[tolIP]]" 
				<< "  [[contactpqrfname]] " << endl; exit(0);
	}
	
	const char* pqrfname   = argv[1];			// arg 1 = atompqrfile ( in )
	const char* vertfname = argv[2];			// arg 2 = vertfile ( in )
	const char* cenfname = argv[3];				// arg 3 = centerfile ( out )
	const double tolSP = atof( argv[4] ); // arg 4 = global tolerance 
																				// (how much sphere can protrude from SES) 
	const int makecontacts = atoi( argv[5] ); // arg 5 = whether contacts are needed
																						// 0 for NO, 1 for YES	
	ofstream fout;
	int np;
	
	// Declaring some variables
	vector<CPnt> scen, apos, capos;
	vector<double> srad, CHG, arad, dum, carad;
	
	// load PQR files
	readpqrAtoms( pqrfname, apos, CHG ,arad );
	
	cout <<"Atoms: "<<apos.size(  )<<endl;  // all atom positions (lab frame)
	
	// read molecular surface points from MSMS vertex files
	vector<CPnt> SP, NP;
	readSurface( vertfname, SP, NP );
	
	// For the contacts make
	const double tolIP = 0.0;
	vector<int> IPind;
	
	// if we need to make contact lists
	if ( makecontacts == 1 ) 
	{
		const double tolIP = atof(argv[6]); // arg 6 = interface tolerance 
		// (same as tolSP, but for interface)
		const char* contactpqrfname = argv[7]; // arg 7 = contact.atom.pqrfile (in)
		
		// load PQR files
		readpqrAtoms( contactpqrfname, capos, dum, carad  );
		
		// find interface contact positions
		vector <CPnt> conPos;
		const double conCutoff = 5.0;
		getInterfaceAtoms( apos, capos, arad, carad, conPos, conCutoff );
		cout <<"Contact Atoms: "<<conPos.size(  )<<endl;
		
		// identify interfacial vertex points
		const double atomVertCutoff = 2.0;
		classifyVertPoints( conPos, SP, IPind, atomVertCutoff  );
	}
	
	// find least number of spheres needed to encompass vertex points
	vector<CPnt> cens;
	vector<double> radii;
	
	if ( makecontacts == 0 ) {findCenters( apos, arad, SP, NP, cens, radii, tolSP );
	} else { findCenters_interface( apos, arad, SP, NP, IPind, cens, radii, tolSP, tolIP );}
	
	/*#########################################################*/
	// extract solvent exposed centers 
	/*#########################################################*/	
	int nk = cens.size(  );
	
	vector<double> rad2( nk );
	for( int ki = 0; ki < nk; ki++ ) 
		rad2[ki] = radii[ki]*radii[ki] *1.0;
	
	vector<int> indSCen;
	for( int ki = 0; ki < nk; ki++ )
	{
		vector<int> neigh;
		findNeighbors( ki, cens, radii, neigh );
		if(  IsSphereExposed(ki, neigh, cens, radii, rad2, SP, NP ) ) 
			indSCen.push_back(ki);
	}
	
	/*#########################################################*/
	// modify centers to ensure they connect
	/*#########################################################*/	
	cout <<"modifying radius to connect"<<endl;
	// Initializing as false
	bool bAllConnected = false;
	const int maxct = 20;
	const REAL incrementR = 0.5;
	int ct = 0;
	
	// While all spheres are not connected and while there have been less than
	// maxct trials	
	while(  !bAllConnected  && ct < maxct  )
	{
		// Set all connected to true
		bAllConnected = true;
		// For each solvent exposed sphere center
		for( int ki=0; ki<indSCen.size( ); ki++ )
		{
			// Print out the sphere
			cout <<"ct "<<ct<<" ki "<<ki<<" -> ";
			// Find its neighbors
			vector<int> neigh;
			findNeighbors( ki, cens, radii, neigh );	
			// If it has none, set bAllConnected to
			// false and increase sphere radius by
			// 0.5 angstroms
			if( neigh.size( ) == 0) 
			{
				bAllConnected = false;
				radii[ki] += incrementR;
				cout <<" incremented R to "<<radii[ki]<<endl;
			}
			// Otherwise it is connected, move on
			else cout <<"connected"<<endl;
		}
		
		ct++;
	}
	
	// Warn if you were unable to connect the spheres
	if( ct==maxct ) 
		cout <<"warning! reached maxct within completely connected all spheres"<<endl;
	
	// Write center j and radius out to file
	fout.open( cenfname );  cout<<"Writing centers to "<<cenfname<<endl;
	for( int ki = 0; ki < indSCen.size( ); ki++)
	{
		int k = indSCen[ki]; 
		fout<<cens[k].x(  )<<" "<<cens[k].y()<<" "
								<<cens[k].z()<<" "<<radii[k]<<endl;
	}
	fout.close(  );
	
	// Calculate the time taken to run whole process
	double end = read_timer(  );
	cout <<"total time taken [s] "<<end-start<<endl;
	
	return 0;   
}	// end 	main1PQR

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
int main( int argc, char ** argv )
{
	return mainMakesphere( argc, argv );
}


