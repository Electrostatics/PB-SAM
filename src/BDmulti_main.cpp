#include "BDmulti.h"

/******************************************************************/
/******************************************************************/
/**
 *  This function runs the core of BDmulti, inputs:
 \param no of thread
 \param salt conc
 \param no. of trajectories
 \param distance for docking contact
 \param number of monomers ( type 0 )
 \param number of monomers ( type 1 )
 \param boxlength
 \param contactfilename
 \param maxtime
 \param runname
 \param system ( 0 = 2s, 1= 107, 2=354, 3=712 )
 \param restartfilename <<<limit to 12
 \param Mol0 pqr file
 \param Mol0 Imat path
 \param Mol0 expansion path
 \param Mol0 expname
 \param Mol1 pqr file
 \param Mol1 Imat path
 \param Mol1 expansion path
 \param Mol1 expname
 ******************************************************************/
int bd_multi( int argc, char ** argv )
{
	int num_threads = atoi( argv[1] );			// Num threads OMP
	const double salt_conc = atof( argv[2] );	// Salt concentration
	const int ntraj = atoi( argv[3] );			// Number of trajectories
	const double consepdist = atof( argv[4] );	// Separation distance [Angstrom]
	const int np1 = atoi( argv[5] );			// Protein 1
	const int np2 = atoi( argv[6] );			// Protein 2
	const double boxlength = atof( argv[7] );	// Length of sim box [Angstr]
	const char* contactfile = argv[8];			// Contact file name
	const double maxtime = atof( argv[9] );		// Max time for simulation
	char* runname = argv[10];					// Name desired for run
	const int systemIndex  = atoi( argv[11] );	// Not realy sure
	
	for( int i=0; i<argc; i++ )
		cout <<argv[i]<<" " << endl;
	cout<<"variables read"<<endl;
	
	// Generating random number seed
	cout.precision( 5 );
	long int idum=( unsigned )time(NULL);
	srand48( idum );
	cout <<"Random Seed = "<<idum<<endl;
	// Defining number of threads for OMP
#ifdef __OMP
	omp_set_num_threads( num_threads );
#endif
	// system parameters NEED UNITSSS
	const double idiel = 4.0;
	const double sdiel = 78.0;
	const double temperature = 298.15;		// Kelvin
	const double kappa = ANGSTROM * sqrt(   (2 * salt_conc * AVOGADRO_NUM / LITRE
																					 * ELECT_CHG * ELECT_CHG )
																			 / ( sdiel* PERMITTIVITY_VAC * KB * temperature  ) );
	
	cout <<"EPS_SOLVENT:"<<sdiel<<" EPS_IN: "<<idiel<<" KAPPA: "<<kappa<<endl;
	// printing out the number of each molecule and the boxlength in angstroms
	cout <<" Np1 ="<<np1<<" Np2 = "<<np2<<" Boxlength [A] ="<<boxlength<<endl;
	// computes the concentration by dividing the # mol1/ boxlen^3
	cout <<" -> Species 0 Paritial Concentration: "<<np1/pow( boxlength,3.0 )
	<<" [A-3]"<<endl;
	// computes the concentration by dividing the # mol2/ boxlen^3
	cout <<" -> Species 1 Paritial Concentration: "<<np2/pow( boxlength,3.0 )
	<<" [A-3]"<<endl;
	
	// Initiating system constants "true" means here we are using periodic boundary
	// condition (S. Liu)
	CSystem::initConstants( kappa, sdiel, temperature, true, boxlength );
	// Initializing the class BDmulti. nmoltype is the number of different molecule types
	// for now it is set at two
	const int nmoltype=2;
	CBDmulti::initConstants( nmoltype );
	// Set separation distance and read in contacts
	CContact::setSepDist( consepdist );
	vector<CMolContact> molcontactlist;
	// Contact information
	int moltype;
	CMolContact::readMolContact( contactfile, moltype, molcontactlist );
	assert ( moltype == 0 ); // since only two species, we always
	// check species 0's contacts
	CBDmulti::addMolContact( moltype, molcontactlist );
	
	// create molecules
	vector<char*> molfnames1, molfnames2;
	molfnames1.resize( 4 );
	molfnames1[0] = argv[12];	//pqr file
	molfnames1[1] = argv[13]; 	//Imat path
	molfnames1[2] = argv[14]; 	//expansion path
	molfnames1[3] = argv[15]; 	//expname
	molfnames2.resize( 4 );
	molfnames2[0] = argv[16];  	//pqr file
	molfnames2[1] = argv[17]; 	//Imat path
	molfnames2[2] = argv[18]; 	//expansion path
	molfnames2[3] = argv[19]; 	//expname
	
	// create BD system
	bool bRestart;
	char* configfname;
	CBDmulti* bd;		// "bd" points to an object of class CBDmulti (S. Liu)
	
	// If there are only 13 arguments, this means that the program
	// is starting from a restart file
	if( argc==13 )
	{
		bRestart = true;
		configfname = argv[12];
		cout <<"Initializing from restart file "<<configfname<<endl;
	}
	// Otherwise, restart from an initial configuration where all
	// molecules are placed on a lattice (S. Liu)
	else
	{
		bRestart = false;
		cout <<"Initializing from lattice "<<endl;
	}
	// create "bd" with a constructor
	bd = new CBDmulti( np1, np2, molfnames1,molfnames2, idiel );
	// start BD
	char outfname[MAX_CHARNUM];
	// Output to go to mfpt.[runname].dat
	sprintf( outfname, "mfpt.%s.dat", runname );
	ofstream fout( outfname );
	// Starting timer
	double start = read_timer(  );
	int cdock = 0;
	
	// For the given number of trajectories
	for ( int r=0; r<ntraj; r++ )
	{
		// Time the run
		double starti = read_timer(  );
		// If a restart configuration is given, use it
		if( bRestart )
			bd->restartConfig( configfname );
		// Otherwise, reset the lattice
		else
			bd->resetLattice(  );		// place molecules on a lattice (S. Liu)
		//  "runFirstDockTime" will perform the BD simulation until
		//  one trajectory ends (S. Liu)
		double docktime = bd->runFirstDockTime(  maxtime, runname  );
		// Time the completed run
		double endi = read_timer(  );
		// If the docktime is less than the max
		// allotted time, add it to the compounds that
		// docked
		if( docktime < maxtime  ) cdock++;
		// print out the simulation time
		fout <<r<<" "<<docktime<<"[ns]; computation time ="<<endi-starti<< endl;
	}
	
	// Report the total docking percent
	cout << "Total Docked "<<cdock<<" out of "<<ntraj<<" trajectories"<<endl;
	// Cleaning up
	CSystem::deleteConstants(  );
	delete bd;
	// Measuring total simulation time
	double end = read_timer(  );
	// Timing reports
	cout <<"Program runtime = "<<end-start<<" s"<< endl;
	cout <<"Program completed successfully"<<endl;
	
	return 0;
} // end BD multi

/******************************************************************/
/******************************************************************/
/**
 * Main for BD multi
 ******************************************************************/
int main( int argc, char ** argv )
{
	return bd_multi( argc, argv );
}

