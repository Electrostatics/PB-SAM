#include "BDnam.h"
#include "omp.h"

/******************************************************************/
/******************************************************************/
/**
 *  	BDNAM portion: initalizes many simulation parameters
 Calls on BDnam.cpp to gather input from external files
 and run BD
 \param no of thread
 \param salt conc
 \param no. of trajectories
 \param distance for docking contact
 \param output filename
 \param runname
 \param contact file location
 \param Mol1 pqr file
 \param Mol1 Imat pat
 \param Mol1 expansion path
 \param Mol1 expname
 \param Mol2 pqr file
 \param Mol2 Imat path
 \param Mol2 expansion path
 \param Mol2 expname
 ******************************************************************/
int bdNam( int argc, char ** argv ) {
	
	// Parameters from input arguments
	int num_threads = atoi( argv[1] );
	const double salt_conc = atof( argv[2] );
	int nTraj = atoi(  argv[3] );
	double consepdist = atof( argv[4] );
	char* outfname = argv[5];
	char* runname = argv[6];
	char* contactfile = argv[7];
	
	// Printing the input arguments //
	for( int i=0; i<argc; i++ )
	{
		cout <<argv[i]<<" ";
		cout<<endl;
	}
	// Getting a random number seed
	cout.precision( 5 );
	// For now set to be the same to monitor simulations
	srand48( 1 );
	cout <<"Random Seed = "<<1<<endl;
	
	// OpenMP parameters, specifically, the number of threads from ARGVs
#ifdef __OMP
	omp_set_num_threads( num_threads );
#endif
	
	// idiel is the dielectric constant for molecule
	const double idiel = 4.0;
	// sdiel is the  dielectric constant for solvent
	const double sdiel = 78.0;
	// the system temperature ( kelvin )
	const double temperature = 298.15;
	// Calculating the inverse debye length, kappa from
	// k = sqrt(  2 * salt conc * Nav / charge1*charge2  ) / ( sdiel * kb * temp )
	const double kappa = ANGSTROM * sqrt(   (2 * salt_conc * AVOGADRO_NUM / LITRE
																					 * ELECT_CHG * ELECT_CHG )
																			 / ( sdiel* PERMITTIVITY_VAC * KB * temperature  ) );
	//Printing out parameters
	cout <<"EPS_SOLVENT:"<<sdiel<<" EPS_IN: "<<idiel<<" KAPPA: "<<kappa<<endl;
	
	// initializing constants for the system, basically sets the kappa and the
	// sdiel and temp in the class CSystem and within that into CMolecule
	// As it is now, there are no periodic BCs
	CSystem::initConstants( kappa, sdiel, temperature ); // non-periodic BC
	
	// Setting nmoltype to 2, not really sure what is going on here
	// something about the contact lists I think
	const int nmoltype=2;
	CBDnam::initConstants( nmoltype );
	
	// read in contacts
	// Sep distance is from ARGV
	CContact::setSepDist( consepdist );
	// Making vector of the contact list
	// List of CG spheres that must be in contact to call
	// the simulation " DOCKED "
	vector<CMolContact> molcontactlist;
	
	// Create an integer moltype
	int moltype;
	// Call object read molcontact
	CMolContact::readMolContact( contactfile, moltype, molcontactlist );
	assert ( moltype == 0 ); // for nam;
	CBDnam::addMolContact( moltype, molcontactlist );
	
	// Creating vectors for the two different molecules
	vector<char*> molfnames1, molfnames2;
	
	// Getting information for MOL1
	molfnames1.resize( 4 );
	molfnames1[0] = argv[8];  //pqr file
	molfnames1[1] = argv[9]; //Imat path
	molfnames1[2] = argv[10]; //expansion path
	molfnames1[3] = argv[11]; //expname
	// Getting information for MOL2
	molfnames2.resize( 4 );
	molfnames2[0] = argv[12];  //pqr file
	molfnames2[1] = argv[13]; //Imat path
	molfnames2[2] = argv[14]; //expansion path
	molfnames2[3] = argv[15]; //expname
	
	CBDnam bd( molfnames1, molfnames2, idiel );
	
	// timing stuff
	ofstream fout( outfname );
	double start = read_timer(  );
	int cdock = 0;
	int scount;
	bool bFirstDock = false;
	
	// Running simulations for ntrajectories
	for ( int i = 0; i < nTraj; i++ )
	{
		double starti = read_timer();
		CBDnam::STATUS status = bd.run( scount, runname );
		
		if ( status == CBDnam::DOCKED )
		{
			cdock++;
			if( !bFirstDock )
			{
				char fname[MAX_CHARNUM];
				sprintf( fname, "%s.mol0.dock.out", runname );
				CMolecule::printMolConfig( bd.getMol( 0 ), fname );
				sprintf( fname, "%s.mol1.dock.out", runname );
				CMolecule::printMolConfig( bd.getMol( 1 ), fname );
				bFirstDock = true;
			}
		}
		// Outputting timing things
		double endi = read_timer();
		fout <<i<<" "<< CBDnam::STATUS_STRING[status] << " " << ( int )status
		<< " ct " << scount <<" traj time ="<<endi-starti<<" s " << endl;
	}
	
	// This is about the number docked per number of trajectories
	// Counts the diffusion or something
	fout <<" "<<cdock << " " << nTraj<<" delta = "<< ( cdock / double(nTraj )) << endl;
	
	// Printing out the total time taken to run the simulation
	double end = read_timer(  );
	printf( "Program run time = %.4lf seconds\n", end-start  );
	return 0;
}	// end bdNam

/******************************************************************/
/******************************************************************/
/**  *  	BDNAM main
 ******************************************************************/
int main( int argc, char ** argv ) {
	
	if ( argc!= 16 )
	{
		cout << "Incorrect usage: ./nam [# threads] [salt conc] ";
		cout << " [# traj] [Dock dist] [ofname] [runname] [contact file]";
		cout << " [Mol1 pqr] [Mol1 Imat] [Mol1 Exp] [Mol1 Expname] ";
		cout << " [Mol2 pqr] [Mol2 Imat] [Mol2 Exp] [Mol2 Expname] "<< endl;
		exit(0);
	}
	return bdNam( argc, argv );
	return 0;
}
