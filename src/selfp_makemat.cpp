#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <stdio.h>
#include <string>
#include <ctime>
#include "molecule.h"
#include "system.h"
#include "readutil.h"

// Initiating some constants,
#define FACT2 5.604586e2													//!< ( COUL_K*IKbT )

const double AVOGADRO_NUM = 6.0221415e23;					//!< Avogadro's number
const double ELECT_CHG = 1.60217646e-19;					//!< The charge of 1 electron in coulombs
const double PERMITTIVITY_VAC = 8.854187817e-12;	//!< the permittivity in a vacuum
const double KB = 1.3806503e-23;									//!< Boltzmann's const
const double ANGSTROM = 1e-10;										//!< Angstrom to meter conversion
const double LITRE = 1e-3;												//!< Litres in kL

using namespace std;                       

/******************************************************************/
/******************************************************************//**
*  main_selfPolarize
	\param flag for self-polarization
	\param an integer for the number of OMP threads
	\param the salt concentration in Moles
	\param the pqr filename for molecule 
	\param imat path a path to the interaction matrices
	\param the name of the selfpolarization run
	\param the interior dielectric (optional)
******************************************************************/
int main_selfPolarize( int argc, char ** argv )
{
	if ( ( argc != 7 ) && ( argc != 8 ) )
	{
		cout << " Correct usage: ./exec spol [nthreads] [Salt conc] [PQR path]" << 
					" [Imat path] [Runname] [[DielP]] " << endl;
	}
	
	cout << "running self-polarization ..." << endl;
	cout.precision( 10 );
	
	// Parallelization
	int num_threads = atoi( argv[2] );
#ifdef __OMP
	omp_set_num_threads( num_threads );
#endif
	
	// Timing start
	double start = read_timer(  );
	
	//////////////// system parameters //////////////
	const double salt_conc = atof( argv[3] );
	const double idiel = (  argc == 7 ? 4.0 :  atof(argv[7] ) );
	const double sdiel = 78;
	const double temp = 298.15;
	const double kappa = ANGSTROM * sqrt( (2 * salt_conc * AVOGADRO_NUM 
																			 / LITRE * ELECT_CHG * ELECT_CHG ) 
																			 / ( sdiel* PERMITTIVITY_VAC * KB * temp ) );
	cout <<"EPS_SOLVENT:"<<sdiel<<" EPS_IN: "<<idiel<<" KAPPA: "<<kappa<<endl;
	
	CSystem::initConstants( kappa, sdiel, temp );
	
	//////////////// read in charge + centers //////////////
	vector<CPnt> scen, POS;
	vector<double> srad, CHG;
	const char * configfile  = argv[4];
	readpqr( configfile, POS, CHG ,srad, scen ); // as PQR 
	
	//////////////// // define center as center of spheres ////////////////
	const int ncen = scen.size(  );// no of unique spheres
	CPnt initialrcen( 0,0,0 );
	for( int i=0; i<ncen; i++ ) 
		initialrcen += scen[i]; initialrcen /= scen.size(  ); 
	cout <<"rcen : " <<initialrcen<<" ; no. of centers:  "<<ncen<<endl;
	
	//////////////// read in interaction matrices //////////////
	vector<REAL*> iMats;
	char* imatpath = argv[5];
	readMats( iMats, N_POLES, imatpath, ncen );
	
	/////////////// generate a copy of spx ///////////////
	vector<vector<CPnt> > SPxes;
	vector<int> nSPx; 
	vector<vector<int> > neighs;
	CMolecule::generateMolSPX( scen, srad, SPxes, nSPx, neighs  );
	
	//////////////// create molecule  //////////////
	const double intraRcutoff = 300;  //set it big so that all spheres are considered near
	
	vector< vector<int> > intraPolLists_far, intraPolLists_near;
	CMolecule::generateMolTypeIntraPolLists( scen, srad, intraRcutoff,
																					intraPolLists_near, intraPolLists_far );
	
	assert ( intraPolLists_near.size( ) == ncen ); // all spheres are 'near' for selfpol
	
	CMolecule mol( initialrcen, scen, srad, CHG, POS, idiel, iMats, 
								intraRcutoff, SPxes, nSPx, neighs, intraPolLists_near );
	double vm, rss;
	process_mem_usage( vm, rss );  cout << "After creating molecules: VM: " 
	<< vm << "; RSS: " << rss << endl;
	
	//////////////// self polarization  //////////////
	bool bPot = true;
	mol.polarize_self( bPot );
	
	// Finish timer, print out total run time
	double end = read_timer(  );
	cout <<"Total Solve Time [s] = "<<end-start<<endl;
	mol.writeMolExpansions( argv[6] );
		
	////////////////  clean up memory    //////////////
	for( int i=0; i<ncen; i++ ) delete [] iMats[i];
	CSystem::deleteConstants(  );
	
	return 0;
} // end main_selfPolarize


/******************************************************************/
/******************************************************************/
/**
*  main_makeImat
	\param flag for iMatrix computation
	\param an integer for the number of OMP threads
	\param the salt concentration in Moles
	\param the pqr filename for molecule 
******************************************************************/
int main_makeImat( int argc, char ** argv )
{
	if (( argc != 5 ) && ( argc != 6 ))
	{
		cout << " Correct usage: ./exec mat [nthreads] [Salt conc] [PQR path] [[DielP]]" << endl;
	}
	cout <<"Building Surface Integral Matrices ..."<<endl;
	
	// Parallel with OMP, get the number of threads
	int num_threads = atoi( argv[2] );
#ifdef __OMP
	omp_set_num_threads( num_threads );
#endif
	
	// Timing run
	double start = read_timer(  );
	
	// system parameters
	const double salt_conc = atof( argv[3] );
	const double idiel = (  argc == 7 ? 4.0 :  atof(argv[5] ) );
	const double sdiel = 78;
	const double temp = 353;
	const double kappa = ANGSTROM * sqrt(   (2 * salt_conc * AVOGADRO_NUM
											 / LITRE * ELECT_CHG * ELECT_CHG )
										 / ( sdiel* PERMITTIVITY_VAC * KB * temp  ) );
	cout <<"EPS_SOLVENT:"<<sdiel<<" EPS_IN: "
	<<idiel<<" KAPPA: "<<kappa<<endl;
	
	CSystem::initConstants( kappa, sdiel, temp );
	
	// read in charge + centers
	vector<CPnt> scen, POS;
	vector<double> srad, CHG;
	const char * configfile  = argv[4];
	readpqr( configfile, POS, CHG ,srad, scen ); // as PQR
	
	// define center as center of spheres
	CPnt rcen( 0,0,0 );
	for( int i=0; i<scen.size( ); i++)
		rcen += scen[i];
	rcen /= scen.size(  );
	cout <<"rcen = "<<rcen<<endl;
	
	/////////////// generate a copy of moltype values ///////////////
	vector<vector<CPnt> > SPxes;
	vector<int> nSPx;
	vector<vector<int> > neighs;
	//// Generate Points for a solvent exposed sphere, only exposed points
	CMolecule::generateMolSPX( scen, srad, SPxes, nSPx, neighs  );
	
	const double intraRcutoff = 300;
	vector< vector<int> > intraPolLists_far, intraPolLists_near;
	///// Generating lists of CG spheres within the mol that are considered near
	// And those that are considered far.
	CMolecule::generateMolTypeIntraPolLists(scen, srad, intraRcutoff,
											intraPolLists_near, intraPolLists_far);
	
	// create molecule
	CMolecule mol(rcen, scen, srad, CHG, POS, idiel, SPxes, nSPx, neighs, intraPolLists_near);
	double end = read_timer(  );
	
	cout <<"Total Time to compute Imat [s] = "<<end-start<<endl;
	
	return 0;
}	// end main_makeImat

/******************************************************************/
/******************************************************************//**
*  MAIN
******************************************************************/
int main( int argc, char ** argv )
{
	// If input arg 1 contains spol, run Selfpolarize
	if ( strncmp(argv[1], "spol", 4 ) == 0)
		return main_selfPolarize( argc, argv );  
	
	// If input arg 1 contains mat, run makeImat
	else if ( strncmp(argv[1], "mat", 3 ) == 0)
		return main_makeImat( argc, argv );  
	
	return 0;
}

