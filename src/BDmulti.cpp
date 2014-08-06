#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <stdio.h>
#include <time.h>

#include "BDmulti.h"
#include "readutil.h"

/*###############################################################################
 * #
 * # File: BDmulti.cpp
 * #
 * # Date: June 2013
 * #
 * # Description: This contains the tools for running
 * # 			BDmulti -  a BD simulation with multiple molecules
 * #
 * # Inputs: 
 * #
 * # Author: EH Yap, L Felberg, S. Liu
 * #
 * # Copyright ( c )
 * #
 * ################################################################################*/

// how often we write restart data to a file
#define WRITEFREQ 500
// the distcutoff_time is the cutoff distance for interactions
#define DISTCUTOFF_TIME 10 // angstroms

using namespace std;

// Listing static variables
REAL CBDmulti::MINDIST;
REAL CBDmulti::MAXDIST;
int CBDmulti::m_nMolType;
vector< vector<CMolContact> > CBDmulti::MOLCONTACTLIST;

/*#########################################################*/
/*#########################################################*/
// static functions:
//
// initConstants
// addMolContact
//
/*#########################################################*/
/*#########################################################*/

/*#########################################################*/
/*#########################################################*/
// A function that initializes the system with cutoffs 
// for interactions and the number of molecules
/*#########################################################*/
/*#########################################################*/

void 
CBDmulti::initConstants(  int nMolType  )
{
	// if the interaction cutoff is smaller that 10A
	if(  DISTCUTOFF_TIME > CMolecule::m_interactRcutoff   )
	{
		// set max dist to 10A
		MAXDIST = DISTCUTOFF_TIME;
		// set min dist to the interaction cutoff 
		MINDIST = CMolecule::m_interactRcutoff;
	} 
	else
	{
		// otherwise, set the max dist to interaction cutoff
		MAXDIST = CMolecule::m_interactRcutoff;  
		// and the min distance to 10A
		MINDIST = DISTCUTOFF_TIME;
	} 

	// Finding the number of different molecules
	m_nMolType = nMolType;
	// for 2 species, we just need one def for the 
	// first molecule contacts
	MOLCONTACTLIST.resize(  1  ); 
} // end initConstants

/*#########################################################*/
/*#########################################################*/
// Creating a contact list
/*#########################################################*/
/*#########################################################*/

void
CBDmulti::addMolContact(  int mol1type, 
				const vector<CMolContact> &molcontactlist  )
{
	MOLCONTACTLIST[0] = molcontactlist ;
}	// end addMolContact

/*#########################################################*/
/*#########################################################*/
// member functions
//
// CBDmulti
/*#########################################################*/
/*#########################################################*/

/*#########################################################*/
/*#########################################################*/
// construct from array placement with random orientation
/*#########################################################*/
/*#########################################################*/

CBDmulti::CBDmulti(  int np1, int np2,    
		const vector<char*> &molfnames1, const vector<char*> &molfnames2, 
		REAL idiel  ) : m_np1(  np1), m_np2(np2)
{
	// Initializing moltypes: 
	initMolTypes(  molfnames1, molfnames2, idiel  );
	
	// Information about process memory usage
	double vm, rss;
	process_mem_usage(  vm, rss  );  
	cout << "VM: " << vm << "; RSS: " << rss << endl;	   

	cout <<"BD initialization complete for "
				<<m_mols.size(    )<<" molecules." <<endl;
}	// end CBDmulti

/*#########################################################*/
/*#########################################################*/
// read in and precompute values for each species
/*#########################################################*/
/*#########################################################*/

void
CBDmulti::initMolTypes(  const vector<char*> &molfnames1,
	  const vector<char*> &molfnames2, REAL idiel  )
{
	// Storing the names of various files for
	// each molecule
	char* pqrfile1  = molfnames1[0];	// PQR
	char* imatpath1 = molfnames1[1];	// Imat
	char* exppath1  = molfnames1[2];	// Self pole path
	char* expname1  = molfnames1[3];	// Selfpole name

	char* pqrfile2  = molfnames2[0];	// PQR
	char* imatpath2 = molfnames2[1];    // Imat
	char* exppath2  = molfnames2[2];    // Self pole path
	char* expname2  = molfnames2[3];    // Selfpole name

	// read in charge + centers
	vector<CPnt> scen1, POS1, scen2, POS2;		// Centers and positions
	vector<double> srad1, CHG1,  srad2, CHG2;	// radius and charges
	readpqr(  pqrfile1, POS1, CHG1 ,srad1, scen1  );
	readpqr(  pqrfile2, POS2, CHG2 ,srad2, scen2  );

	// Something about the size of centers?
	const int ncen1 = scen1.size(    );
	const int ncen2 = scen2.size(    );

	// Setting it to zero intially
	m_initialrcen1 = CPnt(  0,0,0  ), m_initialrcen2 = CPnt(0,0,0);

	// Computing the center of mass of each molecule
	for(  int i=0; i<ncen1; i++  ) 
		m_initialrcen1 += scen1[i]; m_initialrcen1 /= scen1.size(    );
	for(  int i=0; i<ncen2; i++  ) 
		m_initialrcen2 += scen2[i]; m_initialrcen2 /= scen2.size(    );

	// Printing the center of mass
	cout <<"initialrcen1 : " <<m_initialrcen1<<endl;
	cout <<"initialrcen2 : " <<m_initialrcen2<<endl;


	//////////////// read in solved expansions //////////////
	readMats(  m_iMats1, N_POLES, imatpath1, ncen1  );
	readMats(  m_iMats2, N_POLES, imatpath2, ncen2  );
	REAL intraRcutoff_dum;
	readExpansionsP(  m_iF1, m_iH1, N_POLES, exppath1, 
					expname1, ncen1, intraRcutoff_dum  );
	readExpansionsP(  m_iF2, m_iH2, N_POLES, exppath2, 
					expname2, ncen2, intraRcutoff_dum  );

	/////////////// generate values for each moltype ///////////////
	const double intraRcutoff = 30;  //<<<<<<<<
	vector< vector<int> > intraPolLists_far1, intraPolLists_far2;

	// protein 1
	// 'generateMolSPX' will generate surface points ( S. Liu )
	CMolecule::generateMolSPX(  scen1, srad1, m_SPxes1, m_nSPx1, m_neighs1   );
	//'generateMolExposedChargesFH' will calculate surface charge on each sphere
	// represented by multipole expansions of F and H ( S. Liu )
	CMolecule::generateMolExposedChargesFH(  m_iF1, m_iH1, scen1, srad1, m_SPxes1, m_nSPx1,
									 m_qSolvedF1, m_qSolvedH1, m_totalF1, m_totalH1  );
	CMolecule::generateMolCells(  scen1, srad1, m_initialrcen1, m_molcells1   );

	//'generateMolTypeIntraPolLists' will generate lists for 
	//intra-molecular polarization ( S. Liu )
	CMolecule::generateMolTypeIntraPolLists(  scen1, srad1, intraRcutoff,
									  m_intraPolLists_near1, intraPolLists_far1  );

	//see "molecule.cpp" for interpretations of 'computeMolTypeValues'  ( S. Liu )
	CMolecule::computeMolTypeValues(  m_initialrcen1, scen1, srad1, CHG1, POS1, idiel,
							  intraRcutoff,
							  m_SPxes1, m_nSPx1, m_neighs1,
							  m_iF1, m_iH1, m_qSolvedF1, m_qSolvedH1, m_totalF1, m_totalH1,
							  m_intraPolLists_near1,
							  intraPolLists_far1,
							  m_LFs_intraSelf1, m_LHs_intraSelf1,
							  m_LFs_intraSelf_far1, m_LHs_intraSelf_far1  );

	// protein 2
	CMolecule::generateMolSPX(  scen2, srad2, m_SPxes2, m_nSPx2, m_neighs2   );
	CMolecule::generateMolExposedChargesFH(  m_iF2, m_iH2, scen2, srad2, m_SPxes2, m_nSPx2,
									 m_qSolvedF2, m_qSolvedH2, m_totalF2, m_totalH2  );
	CMolecule::generateMolCells(  scen2, srad2, m_initialrcen2, m_molcells2   );

	CMolecule::generateMolTypeIntraPolLists(  scen2, srad2, intraRcutoff,
									  m_intraPolLists_near2, intraPolLists_far2  );

	CMolecule::computeMolTypeValues(  m_initialrcen2, scen2, srad2, CHG2, POS2, idiel,
							  intraRcutoff,
							  m_SPxes2, m_nSPx2, m_neighs2,
							  m_iF2, m_iH2, m_qSolvedF2, m_qSolvedH2, m_totalF2, m_totalH2,
							  m_intraPolLists_near2,
							  intraPolLists_far2,
							  m_LFs_intraSelf2, m_LHs_intraSelf2,
							  m_LFs_intraSelf_far2, m_LHs_intraSelf_far2  );

	// create molecules
	m_mols.clear(    );

	for(  int i=0; i < m_np1; i++  )
	{
		m_mols.push_back(   new CMolecule(0, m_initialrcen1, scen1, srad1, CHG1, POS1, idiel,
					  m_iMats1, intraRcutoff,
					  m_SPxes1, m_nSPx1, m_neighs1, m_intraPolLists_near1,
					  m_iF1, m_iH1, m_LFs_intraSelf1, m_LHs_intraSelf1,
					  m_LFs_intraSelf_far1, m_LHs_intraSelf_far1,
					  m_qSolvedF1, m_qSolvedH1, m_totalF1, m_totalH1, m_molcells1  ));
	}

	for(  int i=0; i < m_np2; i++  )
	{
		m_mols.push_back(   new CMolecule(1, m_initialrcen2, scen2, srad2, CHG2, POS2, idiel,
					  m_iMats2, intraRcutoff,
					  m_SPxes2, m_nSPx2, m_neighs2, m_intraPolLists_near2,
					  m_iF2, m_iH2, m_LFs_intraSelf2, m_LHs_intraSelf2,
					  m_LFs_intraSelf_far2, m_LHs_intraSelf_far2,
					  m_qSolvedF2, m_qSolvedH2, m_totalF2, m_totalH2, m_molcells2  ));
	}

	// Check to make sure the number of spheres is correct
	assert(  m_mols.size(  ) == m_np1+m_np2);

	// Display which molecule indices belong to which moltype
	if(  m_np1 > 0  ) cout << "index 0 to " <<m_np1-1<<" are moltype 0"<<endl;
	if(  m_np2 > 0  ) 
		cout << "index "<<m_np1<<" to "<<m_np1+m_np2-1<<" are moltype 1"<<endl; 

	// =========================================================
	// "interRCutoff", which gives the value of "m_interRCutoff" in class "CMolecule",
	// as you will see in function CMolecule::generateInterXFormsForPolarize_LowMemory,
	// if the distance between sphere i in one molecule and sphere j is less than 
	// 'interRCutoff', we will not polarize the gradients of their expansions, which 
	// will be used for force calculation ( S. Liu )
	const REAL interRCutoff = 10; // <<<<<<<<<<<<<<<<<<<<<<<<
	// interactRCutoff, which gives the valule of "m_interactRCutoff" in class CMolecule
	// is the cutoff radius for doing mutual polarization, i.e., if the distance 
	// between two molecules is greater than interacRCutoff, this pair won't be mutually 
	// polarized ( S. Liu )
	const REAL interactRCutoff = 100; // <<<<<<<<<<<<<<<<<<<<<<<<
	CMolecule::initMutualConstants(  m_mols, interRCutoff, interactRCutoff, true   );
	cout << "interRcutoff = "<<interRCutoff<<endl;
	cout << "interactRcutoff = "<<interactRCutoff<<endl;
	// =========================================================  

	// initialize diffusion parameters
	m_Dtr.resize(  m_nMolType  );
	m_Dr.resize(  m_nMolType  );

	// Diffusion parameters for 1BRS!!!
	m_Dtr[0] = 0.015; m_Dtr[1] = 0.015; // 1BRS values, translational diffusion coefficient
	m_Dr[0] = 4.5e-5; m_Dr[1] = 4e-5;	// rotational diffusion coefficient ( S. Liu )
	return;
} // end initMolTypes

/*#########################################################*/
/*#########################################################*/
// initialize positions, orientation and 
// moltypes on a lattice
/*#########################################################*/
/*#########################################################*/

void 
CBDmulti::resetLattice(    )
{
	// generate all positions and orientation in the array, regardless of moltype
	const int number_of_monomers = m_np1 + m_np2;
	const int number_per_side = int(  ceil( pow(number_of_monomers, 1.0/3.0  )));
	vector<CPnt> rcens;
	vector<CQuat> rots;

	// determine the size of the biggest protein 
	double maxR = DBL_MIN;
	for (  int i=0; i<m_mols.size(  ); i++)
	{ 
		double maxr_i = m_mols[i]->getMaxR(    );
		if(  maxr_i > maxR  ) 
			maxR = maxr_i; 
	}

	// Fact is the lenght given to each molecule in the system
	REAL fact = CSystem::BOXLENGTH / double(  number_per_side  ); 

	// Printing out some simulation details
	cout <<"maxR for all moltype "<<maxR<<endl;
	cout <<"c2c distance between molecules "<<fact<<endl; 
	
	// Checking that no 2 proteins are overlapping
	assert (  fact > 2*maxR  ); 

	// Generating a lattice of positions for the molecules
	const REAL offset = 0.5*CSystem::BOXLENGTH;  
	int n = 0; 
	for(  int i=0; i<number_per_side; i++  )			// X dim
	{
		REAL px = -offset + (  0.5+i  )*fact;
		for(  int j=0; j<number_per_side ; j++  )		// Y dim
		{
			REAL py = -offset + (  0.5+j  )*fact;
			for(  int k=0; k<number_per_side; k++  )	// Z dim
			{
				REAL pz = -offset + (  0.5+k  )*fact;
				rcens.push_back(   CPnt(px,py,pz  ) );
				rots.push_back (    CQuat::chooseRandom(  ) );
			}
		}
	}

	// assign molecules randomly to the generated positions/orientation
	vector<bool> bOccupied(  number_of_monomers, false  );
	// For each molecule
	for(  int i=0; i<number_of_monomers; i++  )
	{
		bool bSpotFound = false;
		// While it doesn't have a spot
		while (   !bSpotFound   )
		{
			// Draw a random number for it
			int k = (  int  ) floor(drand48()*  number_of_monomers );
			// If that position is unoccupied, place it there
			if(   bOccupied[k] == false  ) 
			{
				m_mols[i]->setPos(   rcens[k]   );
				m_mols[i]->setOrient(   rots[k]   ); 
				bSpotFound = true; 
				bOccupied[k] = true;
			}
		}
	}

	// check that nothing collides (  use direct method instead of cell  )
	for(  int i=0; i < number_of_monomers; i++  )
	{
		if(   m_mols[i]->isCollided( m_mols  ) )  
		{
			// Die if molecules collide
			cout <<"Error : initial config clashes"<<i<<endl;
			exit(  1  );
		}
	}
	return;
}

/*#########################################################*/
/*#########################################################*/
//  Program to begin simulation from a given
//  resart file
/*#########################################################*/
/*#########################################################*/

void
CBDmulti::restartConfig(  char* configfname  )
{
	// Allocating space for molecule
	// positions etc
	vector<CPnt> rcens;
	vector<CQuat> rots; 
	vector<int> moltypes, moltypeid, nps;
	readConfigWithMolType(  configfname, rcens, rots, moltypeid, moltypes, nps  ); 

	// for now, two species
	assert (  nps.size(  ) == 2 );
	assert (  nps.size(  ) == moltypes.size() );

	// check that the number of particles agree
	assert (  moltypes[0] == 0 && moltypes[1] == 1  );
	assert (  nps[0] == m_np1 && nps[1] == m_np2  );
	const int number_of_monomers = m_np1 + m_np2;
	assert (  rcens.size(  ) == number_of_monomers && rots.size() == number_of_monomers );

	// set position and orientation  
	for(  int i=0; i<number_of_monomers; i++  )
	{
		m_mols[i]->setPos(   rcens[i]   );
		m_mols[i]->setOrient(   rots[i]   );
	}

	// check that nothing collides (  use direct method instead of cell  )  
	for(  int i=0; i < number_of_monomers; i++  )
	{
		if(   m_mols[i]->isCollided( m_mols  ) )
		{
			// Die if there is a clash in intital configuration
			cout <<"Error : initial config clashes in restart "<<i<<endl;
			exit(  1  );
		}
	}

	return;
}

/*#########################################################*/
/*#########################################################*/
// bd run for 1 trajectory
/*#########################################################*/
/*#########################################################*/

double 
CBDmulti::runFirstDockTime(  double maxtime, char* runname  )
{

	// Initializing a restart file
	char cfgname[200];
	sprintf(  cfgname, "restart_%s.config", runname  );

	// Number of molecules
	int nmol = m_mols.size(    );

	double totaltime = 0; 		// Timer initialization
	double t = 0;				// time in ps
	int n = 0; 					// step number counter
	bool bDocked = false;  		// Status of simulation

	// reset molecule positions and orientations
	// resetLattice(    );

	// Writing out molecules to a file
	CMolecule::writeMolsPQR(  "initial.pqr", m_mols  ); 

	// Starting a timer
	double start = read_timer(    );
	double avgtime = 0.0;				// For the average time of a simulation

	// While the simulation is not Docked and has not exceeded
	// the max allotted time
	while(  !bDocked && t < maxtime   )
	//while(  !bDocked && n < 1  )
	{
		// compute the minimum distance between 2 molecules
		REAL minDist = computeMinSepDist(    );

		// Write out timings if the step number is
		// divisible by the WRITE FREQUENCY
		bool bWrite = (  n % WRITEFREQ == 0  ); 
		if(  bWrite  ) cout <<n<<") Timings ";

		// Timer reading
		double t1 = read_timer(    );

		// Compute the time step size in picoseconds
		// this is the function that does force calculations
		// and moves molecules ( S. Liu )
		REAL dt = propagate(  minDist, bWrite  );  	// in [ps] 
		t += dt;     				// adding it to time

		// for each for the molecules
		for(  int i=0; i<nmol; i++  )
		{
			// Exit loop if this is not molecule of type 0
			if(  m_mols[i]->getMolType(  ) != 0 ) continue;

			// For each of the molecules
			for(  int j=0; j<nmol; j++  )
			{
				// Exit if this is not molecule of type 1
				if(  m_mols[j]->getMolType(  ) != 1 ) continue;
	
				// Check to see if a molecule of type 0 is
				// docked to type 1
				if (  IsDocked(m_mols[i], m_mols[j]  ))
				{
					bDocked = true;
					break;
				}
				// If they are docked, exit this loop, 
				// don't need to check others
				if(  bDocked  ) break;
			}
		}

		double t2 = read_timer(    ); 
		avgtime += t2-t1;

		// If this step neeeds to be written out
		if (  bWrite  )
		{
			cout <<n<<"  ) Timings: dock "<<t2-t1
					<<" Time for this "<<WRITEFREQ<<" steps: "
					<<avgtime/(  int  ) WRITEFREQ<<"s"<< endl;
			avgtime = 0.0;
			cout <<n<< "  ) mindist "<<minDist <<" dt "<<dt<<
			"[ps]; simulated "<<t*PS_TO_NS<<" [ns]"<<endl; 
			CMolecule::saveConfig(  cfgname, m_mols  );
		}
		n++;
	}//end while

	//  CMolecule::writeMolsPQR(  "finished.pqr", m_mols  );

	double end = read_timer(    );
	cout <<"Docked?  " <<bDocked<<" out of "<<n<<
			" steps; comp time for this traj: "<<end-start<<endl;
	
	// Return simulation time in NANOSECONDS
	return t * PS_TO_NS;  
}  // end runFirstDockTime

/*#########################################################*/
/*#########################################################*/
// Function for propagating dynamics w00t
/*#########################################################*/
/*#########################################################*/

double
CBDmulti::propagate(  double minDist, bool bWriteTime  )
{

	// compute minimum distance to decide timestep and whether to compute forces
	bool bComputeForce = (  minDist < CMolecule::m_interactRcutoff   ) ;  
	double dt = compute_dt(  minDist  );

	// compute forces
	const int nmol = m_mols.size(    );
	vector<CPnt> force(  nmol  ), torque( nmol );
	vector<CPnt> dR(  nmol  ), dO( nmol ); 

	#ifndef __DIFF__

		double t1 = read_timer(    );

		if(  bComputeForce  ) 
		{
			// Check to see if compute forces worked
			bool bBlown = CMolecule::computeForces(  m_mols, force, torque  );
			// this is the function that computes forces and torques ( S. Liu )
			
			if(  !bBlown  ) 
			{
				// If it didn't die and report
				cout <<"died somewhere in computeforces"<<endl; 
				CMolecule::saveConfig(  "died.config", m_mols  );
				CMolecule::writeMolsPQR(  "died.pqr", m_mols  );    
				exit(  1  );
			}
		}
	#else
		force.clear(    ); force.resize(nmol);
		torque.clear(    ); torque.resize(nmol);
	#endif

	double t2 = read_timer(    );      
	
	// move molecules
	for(  int i=0; i<nmol; i++  ) 
	{
		int m = m_mols[i]->getMolType(    );
		dR[i] = (  m_Dtr[m]*dt*FACT2  )*force[i];
		dO[i] = (  m_Dr[m]*dt*IKbT*COUL_K  )*torque[i];
	}

	makeMove(  dR, dO, dt  );  		//this is the function that moves molecules ( S. Liu )

	// More timing
	double t3 = read_timer(    );

	// Print out timings to compute forces and move if needed
	if(  bWriteTime  ) 
		cout <<"  computeForces "<<t2-t1<<" move "<<t3-t2<<endl;

	return dt;      
}  // end propagate

/*#########################################################*/
/*#########################################################*/
// Function to move each molecule
/*#########################################################*/
/*#########################################################*/

void
CBDmulti::makeMove(  const vector<CPnt> & dR, const vector<CPnt> & dO, REAL dt  ) 
{

	const int maxc = 500; // for now, stay put if it doesnt move after maxc trials

	// For each molecule
	for(  int i=0; i<m_mols.size(  ); i++)
	{
		int c;	// Counter for number of trials to translate the mol
		// find which type it is
		int m = m_mols[i]->getMolType(    ); 

		// translate
		CPnt dR_; 
		bool bTranslated = false; c=0;			// initiating counter
		while(  !bTranslated && c < maxc  )
		{
			// While the molecule hasn't translated, increment counter
			c++;

			// Compute dR
			// the displacement is the sum of contribution from 
			// electrostatic forces and random displacement, 

			// CBDmulti::getRandVec will generate a random displacement
			// whose distribution is Gaussian with 
			// standard deviation 2*dt*m_Dtr[m] ( S. Liu )
			dR_ = dR[i] + CBDmulti::getRandVec(  sqrt(2*dt*m_Dtr[m]  ));
			// Try to translate
			m_mols[i]->translate(  dR_  );

			// If it is a cell, use this method to check if
			// it collides
			#ifdef __CELL__
				bool bCollided =  m_mols[i]->isCollided_cell(  m_mols  ); 
			#else
				bool bCollided =  m_mols[i]->isCollided(  m_mols  ); 
			#endif
		
			// If it collides, move it back
			if(  bCollided  )	
				m_mols[i]->untranslate(    );
			// Otherwise, keep it moved
			else 
				bTranslated = true;
		}

		// if (  c == maxc  ) { cout << "stuck in translate " << endl; exit(1);}

		// Now rotate 
		CQuat Q;			//CQuat is the class for quarterion ( S. Liu )
		CPnt dO_;
		bool bRotated = false; c = 0;		// setting constants
		// While not rotated
		while(  !bRotated && c < maxc  )
		{
			// Increment counter
			c++;
			// Rotate
			dO_ = dO[i] + CBDmulti::getRandVec(  sqrt(2*dt*m_Dr[m]  ));
			Q = CQuat(  dO_, dO_.norm(  ));

			// Check for collisions
			#ifdef __CELL__
				if(  ! m_mols[i]->willRotCollide_cell(m_mols, Q  ) ) 
			#else
				if(  ! m_mols[i]->willRotCollide(m_mols, Q  ) )
			#endif
			{
				m_mols[i]->rotate(  Q  );	
				bRotated = true;
			}
		}

	// if (  c == maxc  ) 	{ cout << "stuck in rotate " <<i<< endl; exit(1); }
	}
	return; // Ending makeMove
}

/*#########################################################*/
/*#########################################################*/
//  Function to compute the timestep in PS
/*#########################################################*/
/*#########################################################*/

inline REAL 
CBDmulti::compute_dt(  double d  ) const
{
	if (  d - DISTCUTOFF_TIME > 0  )
		return 2.0 + (  d - DISTCUTOFF_TIME  )/15;

	else
		return 2.0;
}

/*#########################################################*/
/*#########################################################*/
// function for calculation of minimum separation 
// between molecuels in a configuration ( S. Liu )
/*#########################################################*/
/*#########################################################*/

// assume DISTCUTOFF_TIME = CMolecule::m_interRcutoff to simplify

REAL 
CBDmulti::computeMinSepDist(    )     
{ 
	double minsepdist = DBL_MAX;  		// The minimum separation dist
	double minM2Msepdist = DBL_MAX;		// The min molecule 2 molecule dist
	const int nmol = m_mols.size(    );	// the number of molecules
	bool bTrackingS2S = false;			// Something about tracking

	// Loop over all pairs of molecules
	for(  int i=0; i<nmol; i++  )
	{
		for(  int j=i+1; j<nmol; j++  )
		{
			// Find their sep distance
			double dM2M = CMolecule::getSepDist(  m_mols[i], m_mols[j]  ); 

			// first check if they are close enough that we care about
			if(  dM2M < MAXDIST   )
			{
				if(  !bTrackingS2S  ) bTrackingS2S = true;

				// evaluate it until it reaches the small limit that we care about
				for(  int ki=0; ki<m_mols[i]->getNKS(  ); ki++)
				{
					for(  int kj=0; kj<m_mols[j]->getNKS(  ); kj++)
					{
						double d = CMolecule::getSepDist(  m_mols[i], ki, m_mols[j], kj  );
						if(   d < minsepdist   )
						{
							if(   d < MINDIST   ) { return MINDIST; }
							else minsepdist = d;
						}
					}// end kj 
				}// end ki
			}

			// otherwise keep track of mol-to-mol distance just in case
			else if (  !bTrackingS2S  )
			{
				if(  dM2M < minM2Msepdist  ) minM2Msepdist  = dM2M;
			}
		}
	} //end (  i,j  )
	

	if(  bTrackingS2S  ) 
		return  minsepdist;
	else  
		return (  minM2Msepdist > 0 ?  minM2Msepdist : 0  ); 

}	// end computeMinSepDist

/*#########################################################*/
/*#########################################################*/
// Checking to see if a protein is docked 
/*#########################################################*/
/*#########################################################*/

// general, for proteins to bind multiple partners                                        
bool
CBDmulti::IsDocked(  CMolecule* mol1, CMolecule* mol2  ) const
{
	// Determine the type of mol1
	int m1 = mol1->getMolType(    );

	// Look up docking definitions for mol1
	int nDockDef = MOLCONTACTLIST[m1].size(    );

	// For each sphere in DockDef
	for(  int h=0; h<nDockDef; h++  )
	{
		CMolContact mcon = MOLCONTACTLIST[m1][h];

		if(  mcon.getMol2Type(  )!= mol2->getMolType() )  continue;
		// assert(  mcon.getMol2Type(  )== mol2->getMolType() ) ;                          

		int ncon = 0;
		vector<CContact> clist = mcon.getContactList(    );

		for(  int k=0; k<clist.size(  ) && ncon < mcon.getNContact(); k++)
		{
			if(   CMolecule::getSepDist(mol1, clist[k].getID1(  ), mol2, clist[k].getID2())
				<= clist[k].getDist(    ) ) ncon++;

			if(   ncon == mcon.getNContact(  )) return true;
		}

	} // end all definitions
	return false;
}

/*#########################################################*/
/*#########################################################*/
// COMMENTED OUT IN OLD CODE!!
/*#########################################################*/
/*#########################################################*/

  /* 
  // debug
  CPnt cen0 = CPnt(  -118.136, 4.67571, 93.6392  ); //i
  CQuat Q0 = CQuat(  -0.679896,CPnt(-0.574849, 0.0429259, -0.453263  ));
  CPnt cen1 = CPnt(  -94.3519, 32.3896, 114.463  ); //j
  CQuat Q1 = CQuat(  -0.11112, CPnt(-0.77886, 0.565594, -0.24725  ));
  m_mols[0]->setPos(  cen0  );
  m_mols[0]->setOrient(  Q0  );
  m_mols[1]->setPos(  cen1  );
  m_mols[1]->setOrient(  Q1  );

  m_mols[0]->rotateRotCoeff(    );
  m_mols[1]->rotateRotCoeff(    );
  
  for(  int ki=0; ki<m_mols[0]->getNKS(  ); ki++)
    for(  int kj=0; kj<m_mols[1]->getNKS(  ); kj++)
      {
	//	cout <<ki<<" "<<kj<<endl;
	m_mols[1]->getKS(  kj  ).rotateExpansions();
	m_mols[0]->reexpandLSFromSelfVal_debug(  m_mols, 0, ki, 1, kj  );
      }
include
  cout << "j>i done"<<endl;

  for(  int ki=0; ki<m_mols[0]->getNKS(  ); ki++)
    for(  int kj=0; kj<m_mols[1]->getNKS(  ); kj++)
      {
        m_mols[0]->getKS(  ki  ).rotateExpansions();
        m_mols[1]->reexpandLSFromSelfVal_debug(  m_mols, 1, kj, 0, ki  );
      }
  cout <<"i>j done"<<endl;
  

  m_mols[1]->getKS(  7  ).rotateExpansions();                                                                   
  m_mols[0]->reexpandLSFromSelfVal_debug(  m_mols, 0, 7, 1, 7  );                                              
  */


/*
// assume DISTCUTOFF_TIME = CMolecule::m_interRcutoff to simplify
REAL 
CBDmulti::computeMinSepDist(    )     
{ 
  double minsepdist = DBL_MAX;
  double minM2Msepdist = DBL_MAX;
  const int nmol = m_mols.size(    );
  bool bTrackingS2S = false;
  //bool bMinReached = false;
  
  //#pragma omp parallel shared(  bMinReached, bTrackingS2S, minsepdist, minM2Msepdist  )
  {
    for(  int i=0; i<nmol; i++  )
      for(  int j=i+1; j<nmol; j++  )
	{
	  double dM2M = CMolecule::getSepDist(  m_mols[i], m_mols[j]  ); 
	  // first check if they are close enough that we care about
	  if(  dM2M < MAXDIST   )
	    {
	      if(  !bTrackingS2S  ) bTrackingS2S = true;
	      // evaluate it until it reaches the small limit that we care about
	      for(  int ki=0; ki<m_mols[i]->getNKS(  ); ki++)
		{
		  //#pragma omp for
		  for(  int kj=0; kj<m_mols[j]->getNKS(  ); kj++)
		    {
		      double d = CMolecule::getSepDist(  m_mols[i], ki, m_mols[j], kj  );
		      if(   d < minsepdist   )
			{
			  if(   d < MINDIST   ) { return MINDIST; }
			  else minsepdist = d;
			}
		    }// end kj (  parallel for  )

		}// end k
	    }
	  // otherwise keep track of mol-to-mol distance just in case
	  else if (  !bTrackingS2S  )
	    {
	      //#pragma omp master 
	      {
		if(  dM2M < minM2Msepdist  ) minM2Msepdist  = dM2M;
	      }
	    }
	  
	}//end (  i,j  )
    
  }// end parallel
  
  if(  bTrackingS2S  ) return  minsepdist;
  else  return (  minM2Msepdist > 0 ?  minM2Msepdist : 0  ); 
  
}
*/

